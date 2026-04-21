// TranscriptToVideo.cxx
//
// Reads a conversation transcript file with [User] / [Claude] turns and
// renders it to an AVI using vtkFFMPEGWriter at 1280x720, 10 fps.
//
// Transcript format:
//   [USER] Some question or command, possibly spanning
//          multiple continuation lines.
//   [ASSISTANT] A response, also possibly
//               spanning multiple lines.
//   ...
//
// Rules:
//   - A line beginning "[User]" or "[Claude]" starts a new turn.
//   - Any subsequent non-blank line that does NOT begin with "[" is a
//     continuation of the previous turn (joined with a single space).
//   - Blank lines are ignored.
//
// Usage:
//   ./TranscriptToVideo <transcript.txt> <output.avi>
//
// Typewriter speed : 10 chars/sec  →  1 frame/char at 10 fps
// Pause after turn : 0.5 s         →  5 frames
// Word-wrap column : 62

#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkCellArray.h>
#include <vtkFFMPEGWriter.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkWindowToImageFilter.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
static const int    WIDTH         = 1280;
static const int    HEIGHT        = 720;
static const int    FPS           = 10;
static const int    CHARS_PER_FRAME = 5;        // 50 chars/sec at 10 fps
static const int    PAUSE_FRAMES  = 3;          // 0.3 s
static const int    WRAP_COLS     = 60;
static const int    MAX_HISTORY   = 18;         // visible completed lines
static const int    FONT_SIZE     = 18;
static const int    TITLE_SIZE    = 23;
static const int    LINE_HEIGHT   = 28;         // pixels between lines
static const int    MARGIN_LEFT   = 48;
static const int    MARGIN_TOP    = 58;         // below title bar

// Colours (R, G, B in 0–1)
static const double COL_BG[3]     = { 0.04, 0.04, 0.10 };   // near-black navy
static const double COL_USER[3]   = { 1.00, 0.84, 0.00 };   // gold
static const double COL_CLAUDE[3] = { 0.68, 0.85, 0.90 };   // light steel blue
static const double COL_TITLE[3]  = { 1.00, 1.00, 1.00 };   // white
static const double COL_BAR[3]    = { 0.30, 0.55, 0.80 };   // progress bar blue
static const double COL_LABEL[3]  = { 0.65, 0.65, 0.65 };   // dim grey for speaker tag

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------
enum class Speaker { User, Claude, Caption };

struct Turn {
    Speaker     speaker;
    std::string text;   // full unwrapped text
};

// A rendered display line (already word-wrapped, tagged with speaker)
struct DisplayLine {
    Speaker     speaker;
    std::string prefix;  // "[USER]      " or "[ASSISTANT] "
    std::string body;    // wrapped segment
    bool        isFirst;   // first line of this turn (prefix printed)
    bool        isBlank;   // spacer between turns (no text rendered)
    bool        isCaption; // caption update — body holds text (empty = clear)
    DisplayLine() : speaker(Speaker::User), isFirst(false), isBlank(false), isCaption(false) {}
};

// ---------------------------------------------------------------------------
// Replace common multi-byte UTF-8 sequences with ASCII equivalents so that
// VTK's font renderer (which only handles ASCII reliably) does not crash.
// Any remaining non-ASCII byte is replaced with '?'.
// ---------------------------------------------------------------------------
static std::string SanitiseASCII(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    const unsigned char* p = (const unsigned char*)s.data();
    const unsigned char* end = p + s.size();
    while (p < end) {
        unsigned char c = *p;
        if (c < 0x80) {
            out += (char)c;
            ++p;
        } else if (c == 0xE2 && (p + 2) < end) {
            // 3-byte sequence starting 0xE2
            unsigned char b1 = *(p+1), b2 = *(p+2);
            // U+2018 ' U+2019 '  (0xE2 0x80 0x98 / 0x99)
            if (b1 == 0x80 && (b2 == 0x98 || b2 == 0x99)) { out += '\''; p += 3; }
            // U+201C " U+201D "  (0xE2 0x80 0x9C / 0x9D)
            else if (b1 == 0x80 && (b2 == 0x9C || b2 == 0x9D)) { out += '"'; p += 3; }
            // U+2013 en-dash (0xE2 0x80 0x93)
            else if (b1 == 0x80 && b2 == 0x93) { out += '-'; p += 3; }
            // U+2014 em-dash (0xE2 0x80 0x94)
            else if (b1 == 0x80 && b2 == 0x94) { out += '-'; p += 3; }
            // U+2026 ellipsis (0xE2 0x80 0xA6)
            else if (b1 == 0x80 && b2 == 0xA6) { out += "..."; p += 3; }
            // U+2192 arrow (0xE2 0x86 0x92)
            else if (b1 == 0x86 && b2 == 0x92) { out += "->"; p += 3; }
            // U+00D7 multiplication sign (0xC3 0x97) — caught below, but
            // U+22xx math operators: replace with '~'
            else { out += '~'; p += 3; }
        } else if (c == 0xC3 && (p + 1) < end) {
            // 2-byte latin-1 supplement
            unsigned char b1 = *(p+1);
            if (b1 == 0x97) { out += 'x'; p += 2; }       // U+00D7 ×
            else if (b1 == 0xB7) { out += '/'; p += 2; }  // U+00F7 ÷
            else { out += '?'; p += 2; }
        } else {
            // Unknown multi-byte: skip the whole sequence
            out += '?';
            ++p;
            while (p < end && (*p & 0xC0) == 0x80) ++p; // skip continuation bytes
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Strip markdown emphasis markers (* and `) so they don't appear in the video.
// ---------------------------------------------------------------------------
static std::string StripMarkdown(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == '*' || c == '`') continue;
        out += c;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Word-wrap: splits `text` into lines of at most maxCols characters.
// ---------------------------------------------------------------------------
static std::vector<std::string> WrapText(const std::string& text, int maxCols)
{
    std::vector<std::string> lines;
    std::istringstream words(text);
    std::string word, line;
    while (words >> word) {
        if (!line.empty() && (int)(line.size() + 1 + word.size()) > maxCols) {
            lines.push_back(line);
            line = word;
        } else {
            if (!line.empty()) line += ' ';
            line += word;
        }
    }
    if (!line.empty()) lines.push_back(line);
    if (lines.empty()) lines.push_back("");
    return lines;
}

// ---------------------------------------------------------------------------
// Parse transcript file into list of Turns.
// ---------------------------------------------------------------------------
static std::vector<Turn> ParseTranscript(const std::string& path)
{
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Cannot open transcript: " << path << "\n";
        std::exit(1);
    }

    std::vector<Turn> turns;
    std::string line;
    while (std::getline(f, line)) {
        // Trim trailing CR (Windows files)
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.empty()) {
            // Blank line within a dialogue turn = paragraph break.
            // Caption turns don't support paragraph breaks.
            if (!turns.empty() && turns.back().speaker != Speaker::Caption) {
                turns.back().text += '\n';   // '\n' sentinel: expands to blank spacer
            }
            continue;
        }

        if (line.size() >= 6 && line.substr(0, 6) == "[USER]") {
            Turn t;
            t.speaker = Speaker::User;
            t.text    = line.substr(6);
            while (!t.text.empty() && t.text.front() == ' ') t.text.erase(t.text.begin());
            t.text = SanitiseASCII(t.text);
            turns.push_back(std::move(t));
        } else if (line.size() >= 11 && line.substr(0, 11) == "[ASSISTANT]") {
            Turn t;
            t.speaker = Speaker::Claude;
            t.text    = line.substr(11);
            while (!t.text.empty() && t.text.front() == ' ') t.text.erase(t.text.begin());
            t.text = SanitiseASCII(t.text);
            turns.push_back(std::move(t));
        } else if (line.size() >= 9 && line.substr(0, 9) == "[CAPTION]") {
            Turn t;
            t.speaker = Speaker::Caption;
            t.text    = line.substr(9);
            while (!t.text.empty() && t.text.front() == ' ') t.text.erase(t.text.begin());
            t.text = SanitiseASCII(t.text);
            turns.push_back(std::move(t));
        } else if (!turns.empty()) {
            // Continuation line — strip leading whitespace, preserve as separate line
            size_t start = line.find_first_not_of(" \t");
            if (start != std::string::npos) {
                if (!turns.back().text.empty()) turns.back().text += '\n';
                turns.back().text += SanitiseASCII(line.substr(start));
            }
        }
    }
    return turns;
}

// ---------------------------------------------------------------------------
// Expand all Turns into DisplayLines (word-wrapped).
// ---------------------------------------------------------------------------
static std::vector<DisplayLine> ExpandToDisplayLines(const std::vector<Turn>& turns)
{
    std::vector<DisplayLine> out;
    // Prefix widths kept equal for alignment
    const std::string userPfx   = "[USER]      ";
    const std::string claudePfx = "[ASSISTANT] ";

    // Track previous non-caption speaker to know when to insert blank spacers
    Speaker prevSpeaker = Speaker::Caption; // sentinel: no spacer before first real turn
    for (size_t ti = 0; ti < turns.size(); ++ti) {
        const Turn& t = turns[ti];

        // Caption turn: emit one special DisplayLine, no spacer, no history entry
        if (t.speaker == Speaker::Caption) {
            DisplayLine dl;
            dl.isCaption = true;
            dl.speaker   = Speaker::Caption;
            // body: whitespace-only → clear; otherwise wrap lines joined with \n
            std::string trimmed = t.text;
            size_t first = trimmed.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) {
                dl.body = "";   // signal to clear caption
            } else {
                // Wrap to ~90 cols for the wider caption zone
                auto lines = WrapText(trimmed, 90);
                for (size_t i = 0; i < lines.size(); ++i) {
                    if (i > 0) dl.body += "\n";
                    dl.body += lines[i];
                }
            }
            out.push_back(dl);
            continue;
        }

        // Insert a blank spacer before each new dialogue turn (not before the first)
        if (prevSpeaker != Speaker::Caption) {
            DisplayLine blank;
            blank.isBlank = true;
            blank.speaker = t.speaker;
            out.push_back(blank);
        }
        prevSpeaker = t.speaker;

        const std::string& pfx = (t.speaker == Speaker::User) ? userPfx : claudePfx;
        int bodyWidth = WRAP_COLS - (int)pfx.size();
        if (bodyWidth < 20) bodyWidth = 20;

        // Split the turn text on '\n' to get individual source lines.
        // Empty source lines become blank spacers (paragraph breaks within a turn).
        // Non-empty source lines are stripped of markdown markers and word-wrapped.
        std::vector<std::string> srcLines;
        {
            std::string seg;
            for (char ch : t.text) {
                if (ch == '\n') { srcLines.push_back(seg); seg.clear(); }
                else seg += ch;
            }
            srcLines.push_back(seg);
        }

        bool firstOfTurn = true;
        for (const auto& srcLine : srcLines) {
            // Strip markdown emphasis markers
            std::string stripped = StripMarkdown(srcLine);
            // Trim leading/trailing whitespace
            size_t f = stripped.find_first_not_of(" \t");
            if (f == std::string::npos) {
                // Blank — paragraph break spacer within the turn
                DisplayLine blank;
                blank.isBlank = true;
                blank.speaker = t.speaker;
                out.push_back(blank);
                continue;
            }
            size_t l = stripped.find_last_not_of(" \t");
            stripped = stripped.substr(f, l - f + 1);

            auto wrapped = WrapText(stripped, bodyWidth);
            for (size_t i = 0; i < wrapped.size(); ++i) {
                DisplayLine dl;
                dl.speaker = t.speaker;
                dl.prefix  = pfx;
                dl.body    = wrapped[i];
                dl.isFirst = firstOfTurn && (i == 0);
                out.push_back(dl);
            }
            firstOfTurn = false;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Helper: set actor text to `s`, configure font/colour, position.
// ---------------------------------------------------------------------------
static void ConfigureTextActor(vtkTextActor* a, int fontSize,
                               double r, double g, double b)
{
    a->GetTextProperty()->SetFontFamilyToCourier();
    a->GetTextProperty()->SetFontSize(fontSize);
    a->GetTextProperty()->SetColor(r, g, b);
    a->GetTextProperty()->SetBold(0);
    a->GetTextProperty()->SetItalic(0);
    a->GetTextProperty()->ShadowOff();
    a->GetTextProperty()->SetJustificationToLeft();
    a->GetTextProperty()->SetVerticalJustificationToBottom();
}

// ---------------------------------------------------------------------------
// Create a thin horizontal progress bar actor (2D polygon).
// Returns the actor; caller adds to renderer.
// ---------------------------------------------------------------------------
static vtkActor2D* MakeProgressBar(vtkPolyData*& pd, vtkPoints*& pts)
{
    pts = vtkPoints::New();
    vtkCellArray* cells = vtkCellArray::New();

    // 4 points, updated every frame
    pts->SetNumberOfPoints(4);
    pts->SetPoint(0, 0,   HEIGHT - 8, 0);
    pts->SetPoint(1, 0,   HEIGHT - 2, 0);
    pts->SetPoint(2, 100, HEIGHT - 2, 0);
    pts->SetPoint(3, 100, HEIGHT - 8, 0);

    vtkIdType quad[4] = { 0, 1, 2, 3 };
    cells->InsertNextCell(4, quad);

    pd = vtkPolyData::New();
    pd->SetPoints(pts);
    pd->SetPolys(cells);
    cells->Delete();

    vtkPolyDataMapper2D* mapper = vtkPolyDataMapper2D::New();
    mapper->SetInputData(pd);

    vtkActor2D* actor = vtkActor2D::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(COL_BAR[0], COL_BAR[1], COL_BAR[2]);
    mapper->Delete();
    return actor;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: TranscriptToVideo <transcript.txt> <output.avi>\n";
        return 1;
    }
    const std::string transcriptPath = argv[1];
    const std::string outputPath     = argv[2];

    // -----------------------------------------------------------------------
    // Parse & expand
    // -----------------------------------------------------------------------
    auto turns        = ParseTranscript(transcriptPath);
    auto displayLines = ExpandToDisplayLines(turns);
    const int totalLines = (int)displayLines.size();

    if (totalLines == 0) {
        std::cerr << "Transcript has no parseable turns.\n";
        return 1;
    }

    // Total frames (rough estimate for progress bar denominator)
    // Each display line: on average ~40 chars + PAUSE_FRAMES
    // We'll compute exact total below after knowing char counts.
    int totalFrames = 0;
    for (const auto& dl : displayLines) {
        int chars = (int)dl.body.size();
        totalFrames += (chars + CHARS_PER_FRAME - 1) / CHARS_PER_FRAME + 1 + PAUSE_FRAMES;
    }

    // -----------------------------------------------------------------------
    // VTK Setup
    // -----------------------------------------------------------------------
    vtkNew<vtkRenderer>     ren;
    vtkNew<vtkRenderWindow> renWin;

    ren->SetBackground(COL_BG[0], COL_BG[1], COL_BG[2]);
    renWin->AddRenderer(ren);
    renWin->SetSize(WIDTH, HEIGHT);
    renWin->SetOffScreenRendering(1);
    renWin->SetMultiSamples(0);

    // -----------------------------------------------------------------------
    // Title actor
    // -----------------------------------------------------------------------
    vtkNew<vtkTextActor> titleActor;
    titleActor->SetInput("Transcript");
    ConfigureTextActor(titleActor, TITLE_SIZE,
                       COL_TITLE[0], COL_TITLE[1], COL_TITLE[2]);
    titleActor->GetTextProperty()->SetBold(1);
    titleActor->GetTextProperty()->SetJustificationToCentered();
    titleActor->SetDisplayPosition(WIDTH / 2, HEIGHT - 36);
    ren->AddActor2D(titleActor);

    // Thin separator line below title (drawn as a narrow quad)
    {
        vtkPoints*    lpts  = vtkPoints::New();
        vtkCellArray* lcell = vtkCellArray::New();
        lpts->InsertNextPoint(0,     HEIGHT - 44, 0);
        lpts->InsertNextPoint(WIDTH, HEIGHT - 44, 0);
        lpts->InsertNextPoint(WIDTH, HEIGHT - 42, 0);
        lpts->InsertNextPoint(0,     HEIGHT - 42, 0);
        vtkIdType q[4] = {0,1,2,3};
        lcell->InsertNextCell(4, q);
        vtkPolyData* lpd = vtkPolyData::New();
        lpd->SetPoints(lpts);
        lpd->SetPolys(lcell);
        vtkPolyDataMapper2D* lmap = vtkPolyDataMapper2D::New();
        lmap->SetInputData(lpd);
        vtkActor2D* lineActor = vtkActor2D::New();
        lineActor->SetMapper(lmap);
        lineActor->GetProperty()->SetColor(0.3, 0.3, 0.5);
        ren->AddActor2D(lineActor);
        lpts->Delete(); lcell->Delete(); lpd->Delete();
        lmap->Delete(); lineActor->Delete();
    }

    // -----------------------------------------------------------------------
    // History text actors (MAX_HISTORY lines, reused with shifting content)
    // -----------------------------------------------------------------------
    std::vector<vtkTextActor*> histActors(MAX_HISTORY);
    for (int i = 0; i < MAX_HISTORY; ++i) {
        histActors[i] = vtkTextActor::New();
        histActors[i]->SetInput("");
        ConfigureTextActor(histActors[i], FONT_SIZE, 0.5, 0.5, 0.5);
        // Position: bottom of history area = just above current-line area
        int y = MARGIN_TOP + (MAX_HISTORY - 1 - i) * LINE_HEIGHT;
        histActors[i]->SetDisplayPosition(MARGIN_LEFT, HEIGHT - y);
        ren->AddActor2D(histActors[i]);
    }

    // Current (actively typing) line actor — two actors: prefix + body
    vtkNew<vtkTextActor> curPfxActor;   // "[User]  " in dim colour
    vtkNew<vtkTextActor> curBodyActor;  // body text in vivid colour
    vtkNew<vtkTextActor> cursorActor;   // blinking block

    ConfigureTextActor(curPfxActor,  FONT_SIZE, COL_LABEL[0], COL_LABEL[1], COL_LABEL[2]);
    ConfigureTextActor(curBodyActor, FONT_SIZE, 1, 1, 1);  // will be overridden per turn
    ConfigureTextActor(cursorActor,  FONT_SIZE, 1, 1, 1);

    cursorActor->SetInput("|");

    int curLineY = HEIGHT - (MARGIN_TOP + MAX_HISTORY * LINE_HEIGHT);
    // Place current line below history
    curPfxActor->SetDisplayPosition(MARGIN_LEFT, curLineY);
    ren->AddActor2D(curPfxActor);
    ren->AddActor2D(curBodyActor);
    ren->AddActor2D(cursorActor);

    // -----------------------------------------------------------------------
    // Progress bar
    // -----------------------------------------------------------------------
    vtkPolyData* barPD  = nullptr;
    vtkPoints*   barPts = nullptr;
    vtkActor2D*  barActor = MakeProgressBar(barPD, barPts);
    ren->AddActor2D(barActor);

    // Progress label
    vtkNew<vtkTextActor> progLabelActor;
    progLabelActor->SetInput("");
    ConfigureTextActor(progLabelActor, 13,
                       COL_LABEL[0], COL_LABEL[1], COL_LABEL[2]);
    progLabelActor->SetDisplayPosition(MARGIN_LEFT, 12);
    ren->AddActor2D(progLabelActor);

    // -----------------------------------------------------------------------
    // Caption actor — bottom-centre, Arial, white, initially hidden
    // -----------------------------------------------------------------------
    vtkNew<vtkTextActor> captionActor;
    captionActor->SetInput("");
    captionActor->GetTextProperty()->SetFontFamilyToArial();
    captionActor->GetTextProperty()->SetFontSize(20);
    captionActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    captionActor->GetTextProperty()->SetBold(0);
    captionActor->GetTextProperty()->SetItalic(1);
    captionActor->GetTextProperty()->ShadowOn();
    captionActor->GetTextProperty()->SetShadowOffset(1, -1);
    captionActor->GetTextProperty()->SetJustificationToCentered();
    captionActor->GetTextProperty()->SetVerticalJustificationToBottom();
    // Position: centred horizontally, in the gap between typing line and progress label
    captionActor->SetDisplayPosition(WIDTH / 2, 32);
    ren->AddActor2D(captionActor);

    // -----------------------------------------------------------------------
    // FFMPEG writer
    // -----------------------------------------------------------------------
    vtkNew<vtkWindowToImageFilter> w2i;
    w2i->SetInput(renWin);
    w2i->SetScale(1);
    w2i->ReadFrontBufferOff();

    vtkNew<vtkFFMPEGWriter> writer;
    writer->SetInputConnection(w2i->GetOutputPort());
    writer->SetFileName(outputPath.c_str());
    writer->SetRate(FPS);
    writer->SetBitRate(4000);
    writer->SetBitRateTolerance(400);
    writer->Start();

    // -----------------------------------------------------------------------
    // Helper: measure pixel width of a string in the current font
    // (approximate: Courier is monospace, so width ≈ chars × charWidth)
    // At font size 17 in Courier, one character ≈ 10.2px wide.
    // -----------------------------------------------------------------------
    const double CHAR_PX = 10.2;

    auto bodyX = [&](const std::string& pfx) -> int {
        return MARGIN_LEFT + (int)(pfx.size() * CHAR_PX);
    };

    // -----------------------------------------------------------------------
    // Render one frame
    // -----------------------------------------------------------------------
    int frameCount = 0;
    auto renderFrame = [&]() {
        renWin->Render();
        w2i->Modified();
        writer->Write();
        ++frameCount;
    };

    // -----------------------------------------------------------------------
    // History ring — holds completed display lines (most-recent last)
    // -----------------------------------------------------------------------
    std::vector<DisplayLine> history;

    auto refreshHistory = [&]() {
        // slot 0 = bottom row (newest), slot MAX_HISTORY-1 = top row (oldest).
        // Brightness fades linearly from 0.85 (slot 0) to 0.20 (slot MAX-1).
        // Blank spacer lines show as empty strings.
        // The [USER]/[ASSISTANT] prefix stays at a fixed bright level so the
        // speaker is always identifiable even in dim history.
        int n = (int)history.size();
        for (int slot = 0; slot < MAX_HISTORY; ++slot) {
            int idx = n - 1 - slot;   // slot 0 → newest, slot MAX-1 → oldest
            if (idx < 0) {
                histActors[slot]->SetInput("");
                continue;
            }
            const auto& hl = history[idx];
            if (hl.isBlank) {
                histActors[slot]->SetInput("");
                continue;
            }

            // Graduated brightness: bright near bottom, dim near top
            double bodyBright = 0.20 + 0.65 * (1.0 - (double)slot / (MAX_HISTORY - 1));

            const double* col = (hl.speaker == Speaker::User) ? COL_USER : COL_CLAUDE;

            if (hl.isFirst) {
                // Prefix stays vivid; body fades
                // We render prefix + body as one string but colour the whole line
                // at body brightness — a compromise since VTK TextActor is single-colour.
                // Use a slightly higher floor for the prefix line so the tag is readable.
                double pfxBright = std::min(1.0, bodyBright + 0.25);
                std::string txt = hl.prefix + hl.body;
                histActors[slot]->SetInput(txt.c_str());
                histActors[slot]->GetTextProperty()->SetColor(
                    col[0] * pfxBright, col[1] * pfxBright, col[2] * pfxBright);
            } else {
                std::string txt = std::string(hl.prefix.size(), ' ') + hl.body;
                histActors[slot]->SetInput(txt.c_str());
                histActors[slot]->GetTextProperty()->SetColor(
                    col[0] * bodyBright, col[1] * bodyBright, col[2] * bodyBright);
            }
        }
    };

    // -----------------------------------------------------------------------
    // Main animation loop
    // -----------------------------------------------------------------------
    int dlIdx = 0;  // index into displayLines

    // Count total turns to feed into progress label
    // (display line index → turn index: count isFirst lines)
    int totalTurns = 0;
    for (const auto& dl : displayLines) if (dl.isFirst) ++totalTurns;
    int turnsSeen = 0;

    for (int li = 0; li < totalLines; ++li) {
        const DisplayLine& dl = displayLines[li];

        // Caption update: swap the bottom caption, no history, no typewriter
        if (dl.isCaption) {
            captionActor->SetInput(dl.body.c_str());
            renderFrame();
            continue;
        }

        // Blank spacer: push to history silently with no typewriter frames
        if (dl.isBlank) {
            history.push_back(dl);
            refreshHistory();
            continue;
        }

        // Vivid colour for current speaker
        double cr, cg, cb;
        if (dl.speaker == Speaker::User) {
            cr = COL_USER[0]; cg = COL_USER[1]; cb = COL_USER[2];
        } else {
            cr = COL_CLAUDE[0]; cg = COL_CLAUDE[1]; cb = COL_CLAUDE[2];
        }
        curBodyActor->GetTextProperty()->SetColor(cr, cg, cb);
        cursorActor->GetTextProperty()->SetColor(cr, cg, cb);

        // Prefix: vivid speaker colour for first line, indent for continuation
        if (dl.isFirst) {
            curPfxActor->SetInput(dl.prefix.c_str());
            curPfxActor->GetTextProperty()->SetColor(cr, cg, cb);
        } else {
            curPfxActor->SetInput(std::string(dl.prefix.size(), ' ').c_str());
            curPfxActor->GetTextProperty()->SetColor(0, 0, 0);
        }

        int bx = bodyX(dl.prefix);
        curBodyActor->SetDisplayPosition(bx, curLineY);
        cursorActor->SetDisplayPosition(bx, curLineY);

        if (dl.isFirst) ++turnsSeen;

        // Update progress label
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "Turn %d / %d", turnsSeen, totalTurns);
            progLabelActor->SetInput(buf);
        }

        // Update progress bar width
        {
            double frac = (totalFrames > 0) ? (double)frameCount / totalFrames : 0.0;
            double barW = frac * (WIDTH - 2 * MARGIN_LEFT);
            barPts->SetPoint(0, MARGIN_LEFT,        HEIGHT - 8, 0);
            barPts->SetPoint(1, MARGIN_LEFT,        HEIGHT - 2, 0);
            barPts->SetPoint(2, MARGIN_LEFT + barW, HEIGHT - 2, 0);
            barPts->SetPoint(3, MARGIN_LEFT + barW, HEIGHT - 8, 0);
            barPts->Modified();
            barPD->Modified();
        }

        // Typewriter: reveal CHARS_PER_FRAME characters per rendered frame.
        // Always render the fully-complete state last.
        const std::string& body = dl.body;
        int ci = 0;
        while (true) {
            curBodyActor->SetInput(body.substr(0, ci).c_str());
            double cx = bx + ci * CHAR_PX;
            cursorActor->SetDisplayPosition((int)cx, curLineY);
            cursorActor->SetVisibility((frameCount % 2 == 0) ? 1 : 0);
            renderFrame();
            if (ci >= (int)body.size()) break;
            ci = std::min(ci + CHARS_PER_FRAME, (int)body.size());
        }

        // Line is complete — move it to history
        history.push_back(dl);
        refreshHistory();

        // Clear current line display
        curPfxActor->SetInput("");
        curBodyActor->SetInput("");
        cursorActor->SetVisibility(0);

        // Pause frames
        for (int p = 0; p < PAUSE_FRAMES; ++p) {
            // Update progress bar
            {
                double frac = (totalFrames > 0) ? (double)frameCount / totalFrames : 0.0;
                double barW = frac * (WIDTH - 2 * MARGIN_LEFT);
                barPts->SetPoint(0, MARGIN_LEFT,        HEIGHT - 8, 0);
                barPts->SetPoint(1, MARGIN_LEFT,        HEIGHT - 2, 0);
                barPts->SetPoint(2, MARGIN_LEFT + barW, HEIGHT - 2, 0);
                barPts->SetPoint(3, MARGIN_LEFT + barW, HEIGHT - 8, 0);
                barPts->Modified();
                barPD->Modified();
            }
            renderFrame();
        }
    }

    // -----------------------------------------------------------------------
    // Finish
    // -----------------------------------------------------------------------
    writer->End();

    // Tidy up manual ref-counted objects
    barActor->Delete();
    barPD->Delete();
    barPts->Delete();
    for (auto* a : histActors) a->Delete();

    std::cerr << "Wrote " << frameCount << " frames ("
              << frameCount / FPS << " s) to " << outputPath << "\n";
    return 0;
}
