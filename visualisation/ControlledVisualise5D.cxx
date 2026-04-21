// ControlledVisualise5D.cxx
// Derived from Visualise5D.cxx by Peter Boyle
//
// A minimal-protocol rendering engine for 5D DWF eigenvector-density data.
// Intended to be driven by an external intelligent controller (e.g. Claude)
// that handles all natural-language interpretation and state tracking.
//
// Commands are sent one per line to the named pipe /tmp/visualise_cmd.
// State is reported to stdout after every command.
//
// Wire protocol (all fields whitespace-separated):
//
//   slice <dim> <N>       set Slice[dim] = N  (0-based, wraps to lattice size)
//   slice <dim> +<N>      increment Slice[dim] by N
//   slice <dim> -<N>      decrement Slice[dim] by N
//   zoom <factor>         camera Dolly by factor  (>1 = in, <1 = out)
//   iso <value>           set isosurface threshold to <value> x RMS
//   file <index>          jump to file by absolute index
//   file +<N>             advance N files
//   file -<N>             go back N files
//   render                force a render with current state
//   status                print current state to stdout
//   quit                  exit cleanly
//
// Dimension indices for 5D DWF grid (e.g. --grid 48.32.32.32.32):
//   s=0 (Ls)  x=1  y=2  z=3  t=4
// For a 4D grid (--grid 32.32.32.32):
//   x=0  y=1  z=2  t=3

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkStripper.h>
#include <vtkImageData.h>
#include <vtkCallbackCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkProperty2D.h>
#include <vtkWindowToImageFilter.h>

#define MPEG
#ifdef MPEG
#include <vtkFFMPEGWriter.h>
#endif

#include <array>
#include <string>
#include <vector>
#include <queue>
#include <mutex>
#include <thread>
#include <atomic>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <Grid/Grid.h>

#define USE_FLYING_EDGES
#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
typedef vtkFlyingEdges3D isosurface;
#else
#include <vtkMarchingCubes.h>
typedef vtkMarchingCubes isosurface;
#endif

#define CMD_PIPE "/tmp/visualise_cmd"

static int g_mpeg      = 0;
static int g_framerate = 10;

// ─── Thread-safe command queue ────────────────────────────────────────────────
static std::queue<std::string> g_cmdQueue;
static std::mutex              g_cmdMutex;
static std::atomic<bool>       g_running{true};
static double                  g_spinDeg = 0.0;   // degrees per poll tick; 0 = stopped

// ─── MPEG recording state ─────────────────────────────────────────────────────
static bool                        g_recording   = false;
static vtkFFMPEGWriter*            g_mpegWriter  = nullptr;
static vtkWindowToImageFilter*     g_imageFilter = nullptr;
static std::string                 g_recordingFile;   // AVI filename for mux step

// ─── Audio state (PCM audio track, synced to video frames) ───────────────────
static const int    AUDIO_RATE        = 44100;
static const double BEEP_FREQ         = 800.0;
static const int    BEEP_SAMPLES      = AUDIO_RATE * 4 / 100;  // 40ms beep
static std::vector<int16_t> g_audioBuffer;
static int          g_beepRemaining   = 0;
static double       g_beepPhase       = 0.0;
static int          g_samplesPerFrame = AUDIO_RATE / 10; // updated at record start

// Write one video frame worth of audio samples (beep or silence) to the buffer.
static void GenerateAudioFrame()
{
    for (int i = 0; i < g_samplesPerFrame; i++) {
        int16_t s = 0;
        if (g_beepRemaining > 0) {
            int pos = BEEP_SAMPLES - g_beepRemaining;
            double env = 1.0;
            int fade = AUDIO_RATE / 100;  // 10ms fade
            if (pos   < fade) env = (double)pos   / fade;
            if (g_beepRemaining < fade) env = (double)g_beepRemaining / fade;
            s = (int16_t)(16000.0 * env * std::sin(2.0 * M_PI * BEEP_FREQ * g_beepPhase / AUDIO_RATE));
            g_beepPhase += 1.0;
            --g_beepRemaining;
        } else {
            g_beepPhase = 0.0;
        }
        g_audioBuffer.push_back(s);
    }
}

static void TriggerBeep() { g_beepRemaining = BEEP_SAMPLES; }

// Simple mono 16-bit PCM WAV writer.
static void WriteWAV(const std::string& path, const std::vector<int16_t>& buf, int rate)
{
    std::ofstream f(path, std::ios::binary);
    int dataBytes  = (int)(buf.size() * 2);
    int chunkSize  = 36 + dataBytes;
    int byteRate   = rate * 2;
    f.write("RIFF", 4); f.write((char*)&chunkSize, 4);
    f.write("WAVE", 4);
    f.write("fmt ", 4);
    int fmtSz = 16;       f.write((char*)&fmtSz,  4);
    int16_t pcm = 1;      f.write((char*)&pcm,    2);
    int16_t ch  = 1;      f.write((char*)&ch,     2);
    f.write((char*)&rate,     4);
    f.write((char*)&byteRate, 4);
    int16_t blk = 2;      f.write((char*)&blk, 2);
    int16_t bps = 16;     f.write((char*)&bps, 2);
    f.write("data", 4);   f.write((char*)&dataBytes, 4);
    f.write((char*)buf.data(), dataBytes);
}

// Play a short audible beep on the local machine (non-blocking).
static void PlayBeepAudible()
{
    system("afplay /System/Library/Sounds/Tink.aiff -v 0.4 &");
}

// ─── Grid I/O ─────────────────────────────────────────────────────────────────
template <class T>
void readFile(T& out, const std::string& fname)
{
    Grid::emptyUserRecord record;
    Grid::ScidacReader RD;
    RD.open(fname);
    RD.readScidacFieldRecord(out, record);
    RD.close();
}

using namespace Grid;

// ─── Command reader thread ────────────────────────────────────────────────────
void CommandReaderThread()
{
    mkfifo(CMD_PIPE, 0666);
    std::cout << "[cmd] Listening on " << CMD_PIPE << std::endl;

    while (g_running) {
        int fd = open(CMD_PIPE, O_RDONLY | O_NONBLOCK);
        if (fd < 0) { usleep(200000); continue; }

        int flags = fcntl(fd, F_GETFL);
        fcntl(fd, F_SETFL, flags & ~O_NONBLOCK);

        char buf[4096];
        std::string partial;
        ssize_t n;
        while (g_running && (n = read(fd, buf, sizeof(buf) - 1)) > 0) {
            buf[n] = '\0';
            partial += buf;
            size_t pos;
            while ((pos = partial.find('\n')) != std::string::npos) {
                std::string line = partial.substr(0, pos);
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (!line.empty()) {
                    std::lock_guard<std::mutex> lk(g_cmdMutex);
                    g_cmdQueue.push(line);
                }
                partial = partial.substr(pos + 1);
            }
        }
        close(fd);
    }
}

// ─── FrameUpdater ─────────────────────────────────────────────────────────────
class FrameUpdater : public vtkCallbackCommand
{
public:
    FrameUpdater() : ffile(0), TimerCount(0), old_file(-1), timerId(-2), maxCount(-1) {}

    static FrameUpdater* New() { return new FrameUpdater; }

    int ffile;
    int old_file;
    int timerId;
    int maxCount;

    Coordinate latt;
    Coordinate xyz_dims, xyz_ranges, g_xyz_ranges;
    uint64_t   xyz_vol;
    Coordinate loop_dims, loop_ranges;
    uint64_t   loop_vol;
    Coordinate sum_dims, sum_ranges;
    uint64_t   sum_vol;
    Coordinate slice_dims;
    Coordinate Slice;

    std::vector<std::string> files;
    int        Nd;
    GridBase*  grid;
    Grid::LatticeComplexD* grid_data;

    double rms;

    vtkImageData* imageData    = nullptr;
    vtkTextActor* text         = nullptr;
    isosurface*   posExtractor = nullptr;
    isosurface*   negExtractor = nullptr;

    void SetGrid(GridBase* _grid)
    {
        grid      = _grid;
        Nd        = grid->Nd();
        latt      = grid->GlobalDimensions();
        grid_data = new Grid::LatticeComplexD(grid);
    }

    void SetFiles(std::vector<std::string> list)  { files = list; old_file = -1; }
    void SetSlice(Coordinate _Slice)              { Slice = _Slice; }

    void SetSumDimensions(Coordinate _SumDims)
    {
        sum_dims = _SumDims; sum_ranges = Coordinate(Nd); sum_vol = 1;
        for (int d = 0; d < Nd; d++) { sum_ranges[d] = sum_dims[d] ? latt[d] : 1; sum_vol *= sum_ranges[d]; }
    }

    void SetLoopDimensions(Coordinate _LoopDims)
    {
        loop_dims = _LoopDims; loop_ranges = Coordinate(Nd); loop_vol = 1;
        for (int d = 0; d < Nd; d++) { loop_ranges[d] = loop_dims[d] ? latt[d] : 1; loop_vol *= loop_ranges[d]; }
    }

    void SetDisplayDimensions(Coordinate _xyz_dims)
    {
        xyz_dims = _xyz_dims; g_xyz_ranges = Coordinate(Nd); xyz_ranges = Coordinate(3); xyz_vol = 1;
        for (int d = 0; d < 3; d++) { xyz_ranges[d] = latt[xyz_dims[d]]; xyz_vol *= xyz_ranges[d]; }
        for (int d = 0; d < Nd; d++) {
            g_xyz_ranges[d] = 1;
            for (int dd = 0; dd < 3; dd++) if (xyz_dims[dd] == d) g_xyz_ranges[d] = latt[d];
        }
    }

    void SetSliceDimensions()
    {
        Coordinate sd;
        for (int d = 0; d < Nd; d++) {
            if (g_xyz_ranges[d] > 1 || loop_dims[d] || sum_dims[d]) continue;
            sd.push_back(d);
        }
        slice_dims = sd;
        std::cout << " Slice dimensions: " << slice_dims << std::endl;
    }

    void FillImageData(int loop_idx)
    {
        Coordinate loop_coor;
        Lexicographic::CoorFromIndex(loop_coor, loop_idx, loop_ranges);
        Coordinate xyz_coor(3), g_xyz_coor(Nd), sum_coor(Nd);
        for (uint64_t xyz = 0; xyz < xyz_vol; xyz++) {
            Lexicographic::CoorFromIndex(xyz_coor,   xyz, xyz_ranges);
            Lexicographic::CoorFromIndex(g_xyz_coor, xyz, g_xyz_ranges);
            RealD val = 0.0;
            for (uint64_t si = 0; si < sum_vol; si++) {
                Lexicographic::CoorFromIndex(sum_coor, si, sum_ranges);
                Coordinate site(Nd);
                for (int d = 0; d < Nd; d++)
                    site[d] = (sum_coor[d] + loop_coor[d] + g_xyz_coor[d] + Slice[d]) % latt[d];
                val += real(peekSite(*grid_data, site));
            }
            imageData->SetScalarComponentFromDouble(xyz_coor[0], xyz_coor[1], xyz_coor[2], 0, val);
        }
        imageData->Modified();
    }

    // Reload if needed, fill image, update label, render — no timer advance.
    void ForceRender(vtkRenderWindowInteractor* iren)
    {
        int file = ((TimerCount / (int)loop_vol) + ffile) % (int)files.size();
        if (file != old_file) {
            std::cout << "[render] Loading " << files[file] << std::endl;
            readFile(*grid_data, files[file]);
            old_file = file;
        }
        FillImageData(TimerCount % (int)loop_vol);
        UpdateLabel(file, TimerCount % (int)loop_vol);
        iren->GetRenderWindow()->Render();
    }

    virtual void Execute(vtkObject* caller, unsigned long eventId, void* callData)
    {
        if (vtkCommand::KeyPressEvent == eventId) {
            vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
            std::string key = iren->GetKeySym();
            if (slice_dims.size() > 0) {
                int vert = slice_dims[slice_dims.size() - 1];
                int horz = slice_dims[0];
                if (key == "Up")    Slice[vert] = (Slice[vert] + 1) % latt[vert];
                if (key == "Down")  Slice[vert] = (Slice[vert] + latt[vert] - 1) % latt[vert];
                if (key == "Right") Slice[horz] = (Slice[horz] + 1) % latt[horz];
                if (key == "Left")  Slice[horz] = (Slice[horz] + latt[horz] - 1) % latt[horz];
            }
            if (key == "greater") ffile = (ffile + 1) % (int)files.size();
            if (key == "less")    ffile = (ffile - 1 + (int)files.size()) % (int)files.size();
            ForceRender(iren);
            return;
        }

        if (vtkCommand::TimerEvent == eventId) {
            // timerId == -2: no animation timer (--notime), ignore all timer events
            if (timerId < 0) return;
            int tid = *(reinterpret_cast<int*>(callData));
            if (tid != timerId) return;
            int file = ((TimerCount / (int)loop_vol) + ffile) % (int)files.size();
            if (file != old_file) { readFile(*grid_data, files[file]); old_file = file; }
            FillImageData(TimerCount % (int)loop_vol);
            UpdateLabel(file, TimerCount % (int)loop_vol);
            dynamic_cast<vtkRenderWindowInteractor*>(caller)->GetRenderWindow()->Render();
            ++TimerCount;
            if (TimerCount >= maxCount && timerId > -1)
                dynamic_cast<vtkRenderWindowInteractor*>(caller)->DestroyTimer(timerId);
        }
    }

private:
    int TimerCount;

    void UpdateLabel(int file, int loop_idx)
    {
        Coordinate loop_coor;
        Lexicographic::CoorFromIndex(loop_coor, loop_idx, loop_ranges);
        // Extract tau value from filename (last '_'-delimited field)
        const std::string& path = files[file];
        std::string tau = path.substr(path.rfind('_') + 1);
        std::stringstream ss;
        ss << "tau = " << tau << "\nSlice " << Slice;
        text->SetInput(ss.str().c_str());
    }
};

// ─── Typewriter caption state ─────────────────────────────────────────────────
// User caption (gold, upper line) — cleared on new user: instruction
static std::string g_userCaptionFull;
static size_t      g_userCaptionPos  = 0;
// Claude caption (light blue, lower line) — cleared when user: arrives
static std::string g_claudeCaptionFull;
static size_t      g_claudeCaptionPos = 0;
static int         g_captionTick = 0;
static const int   g_captionRate = 1;   // ticks per character (1 x 100ms = 10 chars/sec)

static std::string WrapText(const std::string& s, int maxCols = 45) {
    std::istringstream words(s);
    std::string word, line, result;
    while (words >> word) {
        if (!line.empty() && (int)(line.size() + 1 + word.size()) > maxCols) {
            result += line + "\n";
            line = word;
        } else {
            if (!line.empty()) line += " ";
            line += word;
        }
    }
    if (!line.empty()) result += line;
    return result;
}

// ─── CommandHandler ───────────────────────────────────────────────────────────
// Minimal parser for the wire protocol. Natural-language interpretation
// is handled externally (by Claude) before commands reach this program.
class CommandHandler : public vtkCallbackCommand
{
public:
    static CommandHandler* New() { return new CommandHandler; }

    FrameUpdater*               fu;
    vtkCamera*                  camera;
    vtkRenderer*                renderer;
    vtkRenderWindowInteractor*  iren;
    vtkTextActor*               captionActor     = nullptr;  // claude (light blue, lower)
    vtkTextActor*               userCaptionActor = nullptr;  // user  (gold,       upper)
    int    pollTimerId     = -1;
    double isosurfaceLevel = 1.0;   // in RMS units

    void CaptureFrame() {
        if (g_recording && g_mpegWriter && g_imageFilter) {
            GenerateAudioFrame();
            g_imageFilter->Modified();
            g_mpegWriter->Write();
        }
    }

    void SetIsosurface(double level)
    {
        isosurfaceLevel = std::max(0.0, std::min(10.0, level));
        fu->posExtractor->SetValue(0,  isosurfaceLevel * fu->rms);
        fu->negExtractor->SetValue(0, -isosurfaceLevel * fu->rms);
        fu->posExtractor->Modified();
        fu->negExtractor->Modified();
    }

    void PrintStatus()
    {
        std::cout << "[status] file       = " << fu->ffile
                  << " : " << fu->files[fu->ffile] << "\n"
                  << "[status] Slice      = " << fu->Slice << "\n"
                  << "[status] latt       = " << fu->latt << "\n"
                  << "[status] isosurface = " << isosurfaceLevel
                  << " x RMS  (" << isosurfaceLevel * fu->rms << ")" << std::endl;
    }

    // Execute one line of the wire protocol.
    void RunLine(const std::string& line)
    {
        std::istringstream iss(line);
        std::string verb;
        if (!(iss >> verb)) return;

        // ── slice <dim> <N|+N|-N> ────────────────────────────────────────────
        if (verb == "slice") {
            int dim; std::string valstr;
            if (!(iss >> dim >> valstr)) { std::cout << "[cmd] slice: expected dim value" << std::endl; return; }
            if (dim < 0 || dim >= fu->Nd) { std::cout << "[cmd] slice: dim out of range" << std::endl; return; }
            int n = (int)fu->latt[dim];
            int newval;
            if (!valstr.empty() && (valstr[0] == '+' || valstr[0] == '-')) {
                int delta = std::stoi(valstr);
                newval = ((fu->Slice[dim] + delta) % n + n) % n;
            } else {
                newval = ((std::stoi(valstr) % n) + n) % n;
            }
            fu->Slice[dim] = newval;
            fu->ForceRender(iren);
            PrintStatus();
        }

        // ── zoom <factor> ────────────────────────────────────────────────────
        else if (verb == "zoom") {
            double factor;
            if (!(iss >> factor)) { std::cout << "[cmd] zoom: expected factor" << std::endl; return; }
            camera->Dolly(factor);
            renderer->ResetCameraClippingRange();
            iren->GetRenderWindow()->Render();
        }

        // ── azimuth <degrees> ────────────────────────────────────────────────
        else if (verb == "azimuth") {
            double deg;
            if (!(iss >> deg)) { std::cout << "[cmd] azimuth: expected degrees" << std::endl; return; }
            camera->Azimuth(deg);
            renderer->ResetCameraClippingRange();
            iren->GetRenderWindow()->Render();
        }

        // ── elevation <degrees> ──────────────────────────────────────────────
        else if (verb == "elevation") {
            double deg;
            if (!(iss >> deg)) { std::cout << "[cmd] elevation: expected degrees" << std::endl; return; }
            camera->Elevation(deg);
            renderer->ResetCameraClippingRange();
            iren->GetRenderWindow()->Render();
        }

        // ── spin <degrees_per_tick> ──────────────────────────────────────────
        // Applies azimuth rotation on every 100ms poll tick. spin 0 stops.
        else if (verb == "spin") {
            double deg;
            if (!(iss >> deg)) { std::cout << "[cmd] spin: expected degrees" << std::endl; return; }
            g_spinDeg = deg;
            std::cout << "[cmd] spin rate = " << g_spinDeg << " deg/tick" << std::endl;
        }

        // ── caption user: <text> / caption claude: <text> / caption ─────────
        // user:   clears both lines, types user text (gold) on upper line.
        // claude: keeps user line, types response (light blue) on lower line.
        // caption alone clears both immediately.
        else if (verb == "caption") {
            std::string rest;
            std::getline(iss, rest);
            if (!rest.empty() && rest[0] == ' ') rest = rest.substr(1);
            if (rest.empty()) {
                g_userCaptionFull = ""; g_userCaptionPos = 0;
                g_claudeCaptionFull = ""; g_claudeCaptionPos = 0;
                g_captionTick = 0;
                if (userCaptionActor) userCaptionActor->SetInput("");
                if (captionActor)     captionActor->SetInput("");
                iren->GetRenderWindow()->Render(); CaptureFrame();
            } else if (rest.substr(0,5) == "user:") {
                // New instruction: clear both, start typing user text
                g_claudeCaptionFull = ""; g_claudeCaptionPos = 0;
                g_userCaptionFull = WrapText(rest); g_userCaptionPos = 0;
                g_captionTick = 0;
                if (userCaptionActor) userCaptionActor->SetInput("");
                if (captionActor)     captionActor->SetInput("");
                iren->GetRenderWindow()->Render(); CaptureFrame();
            } else {
                // claude: or unlabelled — keep user line, type below
                g_claudeCaptionFull = WrapText(rest); g_claudeCaptionPos = 0;
                g_captionTick = 0;
            }
        }

        // ── record start <filename> / record stop ────────────────────────────
        else if (verb == "record") {
#ifdef MPEG
            std::string sub;
            if (!(iss >> sub)) { std::cout << "[cmd] record: expected start <file> or stop" << std::endl; return; }
            if (sub == "stop") {
                if (g_recording && g_mpegWriter) {
                    g_mpegWriter->End();
                    g_mpegWriter->Delete(); g_mpegWriter = nullptr;
                    g_imageFilter->Delete(); g_imageFilter = nullptr;
                    g_recording = false;
                    std::cout << "[cmd] recording stopped: " << g_recordingFile << std::endl;
                    // Write WAV and mux to MP4
                    std::string wavFile = g_recordingFile + ".wav";
                    WriteWAV(wavFile, g_audioBuffer, AUDIO_RATE);
                    g_audioBuffer.clear();
                    std::string mp4File = g_recordingFile;
                    if (mp4File.size() > 4 && mp4File.substr(mp4File.size()-4) == ".avi")
                        mp4File = mp4File.substr(0, mp4File.size()-4) + ".mp4";
                    else
                        mp4File += ".mp4";
                    std::string cmd = "ffmpeg -y -i \"" + g_recordingFile + "\" -i \"" + wavFile +
                                      "\" -c:v copy -c:a aac -shortest \"" + mp4File + "\" 2>/dev/null";
                    int ret = system(cmd.c_str());
                    if (ret == 0) {
                        std::cout << "[cmd] muxed output: " << mp4File << std::endl;
                        unlink(wavFile.c_str());  // clean up intermediate WAV
                    } else {
                        std::cout << "[cmd] mux failed (ffmpeg not found?). WAV kept: " << wavFile << std::endl;
                    }
                } else {
                    std::cout << "[cmd] not recording" << std::endl;
                }
            } else if (sub == "start") {
                std::string fname = "recording.avi";
                iss >> fname;
                if (g_recording) { std::cout << "[cmd] already recording" << std::endl; return; }
                g_recordingFile    = fname;
                g_audioBuffer.clear();
                g_samplesPerFrame  = AUDIO_RATE / std::max(1, g_framerate);
                g_beepRemaining    = 0;
                g_beepPhase        = 0.0;
                g_imageFilter = vtkWindowToImageFilter::New();
                g_imageFilter->SetInput(iren->GetRenderWindow());
                g_imageFilter->SetInputBufferTypeToRGB();
                g_mpegWriter = vtkFFMPEGWriter::New();
                g_mpegWriter->SetFileName(fname.c_str());
                g_mpegWriter->SetRate(g_framerate);
                g_mpegWriter->SetInputConnection(g_imageFilter->GetOutputPort());
                g_mpegWriter->Start();
                g_recording = true;
                std::cout << "[cmd] recording started: " << fname << std::endl;
            } else {
                std::cout << "[cmd] record: unknown subcommand '" << sub << "'" << std::endl;
            }
#else
            std::cout << "[cmd] record: MPEG support not compiled" << std::endl;
#endif
        }

        // ── iso <value> ──────────────────────────────────────────────────────
        else if (verb == "iso") {
            double val;
            if (!(iss >> val)) { std::cout << "[cmd] iso: expected value" << std::endl; return; }
            SetIsosurface(val);
            fu->ForceRender(iren);
            PrintStatus();
        }

        // ── file <index|+N|-N> ───────────────────────────────────────────────
        else if (verb == "file") {
            std::string valstr;
            if (!(iss >> valstr)) { std::cout << "[cmd] file: expected index" << std::endl; return; }
            int n = (int)fu->files.size();
            int newval;
            if (!valstr.empty() && (valstr[0] == '+' || valstr[0] == '-')) {
                int delta = std::stoi(valstr);
                newval = ((fu->ffile + delta) % n + n) % n;
            } else {
                newval = ((std::stoi(valstr) % n) + n) % n;
            }
            fu->ffile    = newval;
            fu->old_file = -1;
            fu->ForceRender(iren);
            PrintStatus();
        }

        // ── render ───────────────────────────────────────────────────────────
        else if (verb == "render") {
            fu->ForceRender(iren);
        }

        // ── status ───────────────────────────────────────────────────────────
        else if (verb == "status") {
            PrintStatus();
        }

        // ── quit ─────────────────────────────────────────────────────────────
        else if (verb == "quit" || verb == "exit") {
            g_running = false;
            iren->TerminateApp();
        }

        else {
            std::cout << "[cmd] Unknown command: '" << line << "'" << std::endl;
        }
    }

    virtual void Execute(vtkObject*, unsigned long eventId, void* callData)
    {
        if (eventId != vtkCommand::TimerEvent) return;
        if (pollTimerId >= 0) {
            int tid = *(reinterpret_cast<int*>(callData));
            if (tid != pollTimerId) return;
        }

        std::vector<std::string> pending;
        {
            std::lock_guard<std::mutex> lk(g_cmdMutex);
            while (!g_cmdQueue.empty()) { pending.push_back(g_cmdQueue.front()); g_cmdQueue.pop(); }
        }
        for (const auto& line : pending) {
            std::cout << "[cmd] >> " << line << std::endl;
            RunLine(line);
            // CaptureFrame() called inside RunLine for caption; for other
            // rendering commands capture here (duplicate Modified() is harmless)
            CaptureFrame();
        }

        // Typewriter: advance one character every g_captionRate ticks.
        // User line types first; claude line starts once user line is complete.
        bool typing = (g_userCaptionPos < g_userCaptionFull.size()) ||
                      (g_claudeCaptionPos < g_claudeCaptionFull.size());
        if (typing) {
            if (++g_captionTick >= g_captionRate) {
                g_captionTick = 0;
                bool rendered = false;
                if (g_userCaptionPos < g_userCaptionFull.size()) {
                    ++g_userCaptionPos;
                    if (userCaptionActor)
                        userCaptionActor->SetInput(g_userCaptionFull.substr(0, g_userCaptionPos).c_str());
                    PlayBeepAudible();
                    TriggerBeep();
                    rendered = true;
                } else if (g_claudeCaptionPos < g_claudeCaptionFull.size()) {
                    ++g_claudeCaptionPos;
                    if (captionActor)
                        captionActor->SetInput(g_claudeCaptionFull.substr(0, g_claudeCaptionPos).c_str());
                    rendered = true;
                }
                if (rendered) {
                    iren->GetRenderWindow()->Render();
                    CaptureFrame();
                }
            }
        }

        // Apply continuous spin (if active) at poll-timer rate
        if (g_spinDeg != 0.0) {
            camera->Azimuth(g_spinDeg);
            renderer->ResetCameraClippingRange();
            iren->GetRenderWindow()->Render();
            CaptureFrame();
        }
    }
};

// ─── main ─────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[])
{
    using namespace Grid;

    Grid_init(&argc, &argv);
    GridLogLayout();

    auto latt_size   = GridDefaultLatt();
    auto simd_layout = GridDefaultSimd(latt_size.size(), vComplex::Nsimd());
    auto mpi_layout  = GridDefaultMpi();
    GridCartesian Grid(latt_size, simd_layout, mpi_layout);

    double default_contour = 1.0;
    std::string arg;

    std::vector<std::string> file_list({"file1","file2","file3","file4",
                                        "file5","file6","file7","file8"});

    if (GridCmdOptionExists(argv, argv+argc, "--files")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--files");
        GridCmdOptionCSL(arg, file_list);
    }
#ifdef MPEG
    if (GridCmdOptionExists(argv, argv+argc, "--mpeg")) g_mpeg = 1;
#endif
    if (GridCmdOptionExists(argv, argv+argc, "--fps")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--fps");
        GridCmdOptionInt(arg, g_framerate);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--isosurface")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--isosurface");
        GridCmdOptionFloat(arg, default_contour);
    }

    int NoTime = 0, Nd = Grid.Nd();
    Coordinate Slice(Nd,0), SumDims(Nd,0), LoopDims(Nd,0), XYZDims({0,1,2});

    if (GridCmdOptionExists(argv, argv+argc, "--slice")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--slice");
        GridCmdOptionIntVector(arg, Slice);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--sum")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--sum");
        GridCmdOptionIntVector(arg, SumDims);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--loop")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--loop");
        GridCmdOptionIntVector(arg, LoopDims);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--xyz")) {
        arg = GridCmdOptionPayload(argv, argv+argc, "--xyz");
        GridCmdOptionIntVector(arg, XYZDims);
    }
    if (GridCmdOptionExists(argv, argv+argc, "--notime")) { NoTime = 1; }

    std::thread cmdThread(CommandReaderThread);
    cmdThread.detach();

    // ── VTK scene ────────────────────────────────────────────────────────────
    vtkNew<vtkNamedColors> colors;
    std::array<unsigned char,4> posColor{{240,184,160,255}}; colors->SetColor("posColor", posColor.data());
    std::array<unsigned char,4> bkg{{51,77,102,255}};        colors->SetColor("BkgColor", bkg.data());

    vtkNew<vtkRenderWindow>           renWin;
    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin);

    int frameCount = (int)file_list.size();
    for (int d = 0; d < Nd; d++) if (LoopDims[d]) frameCount *= latt_size[d];

    vtkNew<vtkCamera> aCamera;
    aCamera->SetViewUp(0,0,-1); aCamera->SetPosition(0,-1000,0); aCamera->SetFocalPoint(0,0,0);
    aCamera->ComputeViewPlaneNormal(); aCamera->Azimuth(30.0); aCamera->Elevation(30.0);

    vtkNew<vtkRenderer> aRenderer;
    renWin->AddRenderer(aRenderer);

    double nrm, rms, contour;
    { LatticeComplexD data(&Grid); readFile(data, file_list[0]); nrm = norm2(data); }
    rms     = std::sqrt(nrm / Grid.gSites());
    contour = default_contour * rms;

    vtkNew<vtkImageData> imageData;
    imageData->SetDimensions(latt_size[XYZDims[0]], latt_size[XYZDims[1]], latt_size[XYZDims[2]]);
    imageData->AllocateScalars(VTK_DOUBLE, 1);
    for (int xx=0;xx<latt_size[XYZDims[0]];xx++)
    for (int yy=0;yy<latt_size[XYZDims[1]];yy++)
    for (int zz=0;zz<latt_size[XYZDims[2]];zz++)
        imageData->SetScalarComponentFromDouble(xx,yy,zz,0,0.0);

    vtkNew<isosurface> posExtractor; posExtractor->SetInputData(imageData); posExtractor->SetValue(0,  contour);
    vtkNew<vtkStripper> posStripper; posStripper->SetInputConnection(posExtractor->GetOutputPort());
    vtkNew<vtkPolyDataMapper> posMapper; posMapper->SetInputConnection(posStripper->GetOutputPort()); posMapper->ScalarVisibilityOff();
    vtkNew<vtkActor> pos; pos->SetMapper(posMapper);
    pos->GetProperty()->SetDiffuseColor(colors->GetColor3d("posColor").GetData());
    pos->GetProperty()->SetSpecular(0.3); pos->GetProperty()->SetSpecularPower(20); pos->GetProperty()->SetOpacity(0.5);

    vtkNew<isosurface> negExtractor; negExtractor->SetInputData(imageData); negExtractor->SetValue(0, -contour);
    vtkNew<vtkStripper> negStripper; negStripper->SetInputConnection(negExtractor->GetOutputPort());
    vtkNew<vtkPolyDataMapper> negMapper; negMapper->SetInputConnection(negStripper->GetOutputPort()); negMapper->ScalarVisibilityOff();
    vtkNew<vtkActor> neg; neg->SetMapper(negMapper);
    neg->GetProperty()->SetDiffuseColor(colors->GetColor3d("Ivory").GetData());

    vtkNew<vtkOutlineFilter> outlineData; outlineData->SetInputData(imageData);
    vtkNew<vtkPolyDataMapper> mapOutline; mapOutline->SetInputConnection(outlineData->GetOutputPort());
    vtkNew<vtkActor> outline; outline->SetMapper(mapOutline);
    outline->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    vtkNew<vtkTextActor> TextT;
    TextT->SetInput("Initialising...");
    TextT->SetPosition(10, 920);
    TextT->GetTextProperty()->SetFontSize(24);
    TextT->GetTextProperty()->SetColor(colors->GetColor3d("Gold").GetData());

    // Claude response caption (light blue, lower line)
    vtkNew<vtkTextActor> CaptionT;
    CaptionT->SetInput("");
    CaptionT->SetPosition(512, 38);
    CaptionT->GetTextProperty()->SetFontSize(32);
    CaptionT->GetTextProperty()->SetColor(0.6, 0.9, 1.0);
    CaptionT->GetTextProperty()->SetJustificationToCentered();
    CaptionT->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
    CaptionT->GetTextProperty()->SetBackgroundOpacity(0.6);
    CaptionT->GetTextProperty()->BoldOn();

    // User instruction caption (gold, upper line)
    vtkNew<vtkTextActor> UserCaptionT;
    UserCaptionT->SetInput("");
    UserCaptionT->SetPosition(512, 82);
    UserCaptionT->GetTextProperty()->SetFontSize(32);
    UserCaptionT->GetTextProperty()->SetColor(1.0, 0.85, 0.0);
    UserCaptionT->GetTextProperty()->SetJustificationToCentered();
    UserCaptionT->GetTextProperty()->SetBackgroundColor(0.0, 0.0, 0.0);
    UserCaptionT->GetTextProperty()->SetBackgroundOpacity(0.6);
    UserCaptionT->GetTextProperty()->BoldOn();

    aRenderer->AddActor(TextT); aRenderer->AddActor(CaptionT); aRenderer->AddActor(UserCaptionT); aRenderer->AddActor(outline);
    aRenderer->AddActor(pos);   aRenderer->AddActor(neg);

    vtkNew<FrameUpdater> fu;
    fu->SetGrid(&Grid); fu->SetFiles(file_list); fu->SetSlice(Slice);
    fu->SetSumDimensions(SumDims); fu->SetLoopDimensions(LoopDims);
    fu->SetDisplayDimensions(XYZDims); fu->SetSliceDimensions();
    fu->imageData = imageData; fu->text = TextT; fu->maxCount = frameCount;
    fu->posExtractor = posExtractor; fu->negExtractor = negExtractor; fu->rms = rms;

    iren->AddObserver(vtkCommand::TimerEvent,    fu);
    iren->AddObserver(vtkCommand::KeyPressEvent, fu);

    aRenderer->SetActiveCamera(aCamera); aRenderer->ResetCamera();
    aRenderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
    aCamera->Dolly(1.0); aRenderer->SetViewport(0.0,0.0,1.0,1.0);
    aRenderer->ResetCameraClippingRange();
    renWin->SetSize(1024,1024); renWin->SetWindowName("ControlledFieldDensity");
    renWin->Render(); iren->Initialize();

    // CommandHandler on fast poll timer
    vtkNew<CommandHandler> cmdHandler;
    cmdHandler->fu             = fu;
    cmdHandler->camera         = aCamera;
    cmdHandler->renderer       = aRenderer;
    cmdHandler->iren           = iren;
    cmdHandler->captionActor     = CaptionT;
    cmdHandler->userCaptionActor = UserCaptionT;
    cmdHandler->isosurfaceLevel = default_contour;
    iren->AddObserver(vtkCommand::TimerEvent, cmdHandler);
    cmdHandler->pollTimerId = iren->CreateRepeatingTimer(100);

    if (g_mpeg == 0 && NoTime == 0) {
        fu->timerId = iren->CreateRepeatingTimer(10000 / g_framerate);
    }

    if (g_mpeg) {
#ifdef MPEG
        vtkWindowToImageFilter* imageFilter = vtkWindowToImageFilter::New();
        imageFilter->SetInput(renWin); imageFilter->SetInputBufferTypeToRGB();
        vtkFFMPEGWriter* writer = vtkFFMPEGWriter::New();
        writer->SetFileName("movie.avi"); writer->SetRate(g_framerate);
        writer->SetInputConnection(imageFilter->GetOutputPort()); writer->Start();
        for (int i = 0; i < fu->maxCount; i++) {
            fu->Execute(iren, vtkCommand::TimerEvent, &fu->timerId);
            imageFilter->Modified(); writer->Write();
        }
        writer->End(); writer->Delete(); imageFilter->Delete();
#else
        assert(-1 && "MPEG support not compiled");
#endif
    } else {
        iren->Start();
    }

    g_running = false;
    Grid_finalize();
    return EXIT_SUCCESS;
}
