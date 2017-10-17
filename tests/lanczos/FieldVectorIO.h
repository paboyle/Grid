namespace Grid { 

 namespace FieldVectorIO {

   // zlib's crc32 gets 0.4 GB/s on KNL single thread
   // below gets 4.8 GB/s
   static uint32_t crc32_threaded(unsigned char* data, int64_t len, uint32_t previousCrc32 = 0) {

     // crc32 of zlib was incorrect for very large sizes, so do it block-wise
     uint32_t crc = previousCrc32;
     off_t blk = 0;
     off_t step = 1024*1024*1024;
     while (len > step) {
       crc = crc32(crc,&data[blk],step);
       blk += step;
       len -= step;
     }

     crc = crc32(crc,&data[blk],len);
     return crc;

   }
 
   static int get_bfm_index( int* pos, int co, int* s ) {
     
     int ls = s[0];
     int NtHalf = s[4] / 2;
     int simd_coor = pos[4] / NtHalf;
     int regu_coor = (pos[1] + s[1] * (pos[2] + s[2] * ( pos[3] + s[3] * (pos[4] % NtHalf) ) )) / 2;
     
     return regu_coor * ls * 48 + pos[0] * 48 + co * 4 + simd_coor * 2;
   }
   
   static void get_read_geometry(const GridBase* _grid,const std::vector<int>& cnodes,
			  std::map<int, std::vector<int> >& slots, 
			  std::vector<int>& slot_lvol,
			  std::vector<int>& lvol,
			  int64_t& slot_lsites,int& ntotal) {

     int _nd = (int)cnodes.size();
     std::vector<int> nodes = cnodes;

     slots.clear();
     slot_lvol.clear();
     lvol.clear();

     int i;
     ntotal = 1;
     int64_t lsites = 1;
     slot_lsites = 1;
     for (i=0;i<_nd;i++) {
       assert(_grid->_fdimensions[i] % nodes[i] == 0);
       slot_lvol.push_back(_grid->_fdimensions[i] / nodes[i]);
       lvol.push_back(_grid->_fdimensions[i] / _grid->_processors[i]);
       lsites *= lvol.back();
       slot_lsites *= slot_lvol.back();
       ntotal *= nodes[i];
     }

     std::vector<int> lcoor, gcoor, scoor;
     lcoor.resize(_nd); gcoor.resize(_nd);  scoor.resize(_nd);
     
     // create mapping of indices to slots
     for (int lidx = 0; lidx < lsites; lidx++) {
       Lexicographic::CoorFromIndex(lcoor,lidx,lvol);
       for (int i=0;i<_nd;i++) {
	 gcoor[i] = lcoor[i] + _grid->_processor_coor[i]*lvol[i];
	 scoor[i] = gcoor[i] / slot_lvol[i];
       }
       int slot;
       Lexicographic::IndexFromCoor(scoor,slot,nodes);
       auto sl = slots.find(slot);
       if (sl == slots.end())
	 slots[slot] = std::vector<int>();
       slots[slot].push_back(lidx);
     }
   }
   
   static void canonical_block_to_coarse_coordinates(GridBase* _coarsegrid,int nb,int& ii,int& oi) {
      // canonical nb needs to be mapped in a coordinate on my coarsegrid (ii,io)
      std::vector<int> _l = _coarsegrid->LocalDimensions();
      std::vector<int> _cl = { _l[1], _l[2], _l[3], _l[4], _l[0] };
      std::vector<int> _cc(_l.size());
      Lexicographic::CoorFromIndex(_cc,nb,_cl);
      std::vector<int> _c = { _cc[4], _cc[0], _cc[1], _cc[2], _cc[3] };
      ii = _coarsegrid->iIndex(_c);
      oi = _coarsegrid->oIndex(_c);
    }

    template<typename Field>
     static bool read_argonne(BasisFieldVector<Field>& ret,const char* dir, const std::vector<int>& cnodes) {

     GridBase* _grid = ret._v[0]._grid;

     std::map<int, std::vector<int> > slots;
     std::vector<int> slot_lvol, lvol;
     int64_t slot_lsites;
     int ntotal;
     get_read_geometry(_grid,cnodes,
		       slots,slot_lvol,lvol,slot_lsites,
		       ntotal);
     int _nd = (int)lvol.size();

     // this is slow code to read the argonne file format for debugging purposes
     int nperdir = ntotal / 32;
     if (nperdir < 1)
       nperdir=1;
     std::cout << GridLogMessage << " Read " << dir << " nodes = " << cnodes << std::endl;
     std::cout << GridLogMessage << " lvol = " << lvol << std::endl;
     
     // for error messages
     char hostname[1024];
     gethostname(hostname, 1024);

     // now load one slot at a time and fill the vector
     for (auto sl=slots.begin();sl!=slots.end();sl++) {
       std::vector<int>& idx = sl->second;
       int slot = sl->first;
       std::vector<float> rdata;
       
       char buf[4096];
       
       sprintf(buf,"%s/checksums.txt",dir); printf("read_argonne: Reading from %s\n",buf);
       FILE* f = fopen(buf,"rt");
       if (!f) {
	 fprintf(stderr,"Node %s cannot read %s\n",hostname,buf); fflush(stderr);
	 return false;
       }
       
       for (int l=0;l<3+slot;l++)
	 fgets(buf,sizeof(buf),f);
       uint32_t crc_exp = strtol(buf, NULL, 16);
       fclose(f);
       
       // load one slot vector
       sprintf(buf,"%s/%2.2d/%10.10d",dir,slot/nperdir,slot);
       f = fopen(buf,"rb");
       if (!f) {
	 fprintf(stderr,"Node %s cannot read %s\n",hostname,buf); fflush(stderr);
	 return false;
       }
       
       fseeko(f,0,SEEK_END);
       off_t total_size = ftello(f);
       fseeko(f,0,SEEK_SET);
       
       int64_t size = slot_lsites / 2 * 24*4;
       rdata.resize(size);
       
       assert(total_size % size == 0);
       
       int _Nfull = total_size / size;
       ret._v.resize(_Nfull,ret._v[0]);
       ret._Nm = _Nfull;
       
       uint32_t crc = 0x0;
       GridStopWatch gsw,gsw2;
       for (int nev = 0;nev < _Nfull;nev++) {
	 
	 gsw.Start();
	 assert(fread(&rdata[0],size,1,f) == 1);
	 gsw.Stop();
	 
	 gsw2.Start();
	 crc = crc32_threaded((unsigned char*)&rdata[0],size,crc);
	 gsw2.Stop();
	 
	 for (int i=0;i<size/4;i++) {
	   char* c = (char*)&rdata[i];
	   char tmp; int j;
	   for (j=0;j<2;j++) {
	     tmp = c[j]; c[j] = c[3-j]; c[3-j] = tmp;
	   }
	 }
	 
	 // loop
	 ret._v[nev].checkerboard = Odd;
#pragma omp parallel 
	 {
	   
	   std::vector<int> lcoor, gcoor, scoor, slcoor;
	   lcoor.resize(_nd); gcoor.resize(_nd); 
	   slcoor.resize(_nd); scoor.resize(_nd);
	   
#pragma omp for
	   for (int64_t lidx = 0; lidx < idx.size(); lidx++) {
	     int llidx = idx[lidx];
	     Lexicographic::CoorFromIndex(lcoor,llidx,lvol);
	     for (int i=0;i<_nd;i++) {
	       gcoor[i] = lcoor[i] + _grid->_processor_coor[i]*lvol[i];
	       scoor[i] = gcoor[i] / slot_lvol[i];
	       slcoor[i] = gcoor[i] - scoor[i]*slot_lvol[i];
	     }
	     
	     if ((lcoor[1]+lcoor[2]+lcoor[3]+lcoor[4]) % 2 == 1) {
	       // poke 
	       iScalar<iVector<iVector<ComplexF, 3>, 4> > sc;
	       for (int s=0;s<4;s++)
		 for (int c=0;c<3;c++)
		   sc()(s)(c) = *(std::complex<float>*)&rdata[get_bfm_index(&slcoor[0],c+s*3, &slot_lvol[0] )];
	       
	       pokeLocalSite(sc,ret._v[nev],lcoor);
	     }
	     
	   }
	 }
       }
       
       fclose(f);      
       std::cout << GridLogMessage << "Loading slot " << slot << " with " << idx.size() << " points and " 
		 << _Nfull << " vectors in "
		 << gsw.Elapsed() << " at " 
		 << ( (double)size * _Nfull / 1024./1024./1024. / gsw.useconds()*1000.*1000. )
		 << " GB/s " << " crc32 = " << std::hex << crc << " crc32_expected = " << crc_exp << std::dec
		 << " computed at "
		 << ( (double)size * _Nfull / 1024./1024./1024. / gsw2.useconds()*1000.*1000. )
		 << " GB/s "
		 << std::endl;
       
       assert(crc == crc_exp);
     }

     _grid->Barrier();
     std::cout << GridLogMessage  << "Loading complete" << std::endl;
     
     return true;
   }

   template<typename Field>
     static bool read_argonne(BasisFieldVector<Field>& ret,const char* dir) {


     GridBase* _grid = ret._v[0]._grid;

     char buf[4096];
     sprintf(buf,"%s/nodes.txt",dir);
     FILE* f = fopen(buf,"rt");
     if (!f) {
       if (_grid->IsBoss()) {
	 fprintf(stderr,"Attempting to load eigenvectors without secifying node layout failed due to absence of nodes.txt\n");
	 fflush(stderr);
       }
       return false;
     }


     std::vector<int> nodes((int)_grid->_processors.size());
     for (int i =0;i<(int)_grid->_processors.size();i++)
       assert(fscanf(f,"%d\n",&nodes[i])==1);
     fclose(f);

     return read_argonne(ret,dir,nodes);
   }

   static void flush_bytes(FILE* f, std::vector<char>& fbuf) {
     if (fbuf.size()) {

       if (fwrite(&fbuf[0],fbuf.size(),1,f) != 1) {
	 fprintf(stderr,"Write failed of %g GB!\n",(double)fbuf.size() / 1024./1024./1024.);
	 exit(2);
       }

       fbuf.resize(0);

     }
   }

   static void write_bytes(void* buf, int64_t s, FILE* f, std::vector<char>& fbuf, uint32_t& crc) {
     static double data_counter = 0.0;
     static GridStopWatch gsw_crc, gsw_flush1,gsw_flush2,gsw_write,gsw_memcpy;
     if (s == 0)
       return;

     // checksum
     gsw_crc.Start();
     crc = crc32_threaded((unsigned char*)buf,s,crc);
     gsw_crc.Stop();

     if (s > fbuf.capacity()) {
       // cannot buffer this, so first flush current buffer contents and then write this directly to file
       gsw_flush1.Start();
       flush_bytes(f,fbuf);
       gsw_flush1.Stop();

       gsw_write.Start();
       if (fwrite(buf,s,1,f) != 1) {
	 fprintf(stderr,"Write failed of %g GB!\n",(double)s / 1024./1024./1024.);
	 exit(2);
       }
       gsw_write.Stop();

     }

     // no room left in buffer, flush to disk
     if (fbuf.size() + s > fbuf.capacity()) {
       gsw_flush2.Start();
       flush_bytes(f,fbuf);
       gsw_flush2.Stop();
     }

     // then fill buffer again
     {
       gsw_memcpy.Start();
       size_t t = fbuf.size();
       fbuf.resize(t + s);
       memcpy(&fbuf[t],buf,s);
       gsw_memcpy.Stop();
     }

     data_counter += (double)s;
     if (data_counter > 1024.*1024.*20.) {
       std::cout << GridLogMessage << "Writing " << ((double)data_counter / 1024./1024./1024.) << " GB at"
	 " crc = " << gsw_crc.Elapsed() << " flush1 = " << gsw_flush1.Elapsed() << " flush2 = " << gsw_flush2.Elapsed() <<
	 " write = " << gsw_write.Elapsed() << " memcpy = " << gsw_memcpy.Elapsed() << std::endl;
       data_counter = 0.0;
       gsw_crc.Reset();
       gsw_write.Reset();
       gsw_memcpy.Reset();
       gsw_flush1.Reset();
       gsw_flush2.Reset();
     }
   }

   static void write_floats(FILE* f, std::vector<char>& fbuf, uint32_t& crc, float* buf, int64_t n) {
     write_bytes(buf,n*sizeof(float),f,fbuf,crc);
   }

   static void read_floats(char* & ptr, float* out, int64_t n) {
     float* in = (float*)ptr;
     ptr += 4*n;

     for (int64_t i=0;i<n;i++)
       out[i] = in[i];
   }

   static int fp_map(float in, float min, float max, int N) {
     // Idea:
     //
     // min=-6
     // max=6
     //
     // N=1
     // [-6,0] -> 0, [0,6] -> 1;  reconstruct 0 -> -3, 1-> 3
     //
     // N=2
     // [-6,-2] -> 0, [-2,2] -> 1, [2,6] -> 2;  reconstruct 0 -> -4, 1->0, 2->4
     int ret =  (int) ( (float)(N+1) * ( (in - min) / (max - min) ) );
     if (ret == N+1) {
       ret = N;
     }
     return ret;
   }

   static float fp_unmap(int val, float min, float max, int N) {
     return min + (float)(val + 0.5) * (max - min)  / (float)( N + 1 );
   }

#define SHRT_UMAX 65535
#define FP16_BASE 1.4142135623730950488
#define FP16_COEF_EXP_SHARE_FLOATS 10
   static float unmap_fp16_exp(unsigned short e) {
     float de = (float)((int)e - SHRT_UMAX / 2);
     return ::pow( FP16_BASE, de );
   }

   // can assume that v >=0 and need to guarantee that unmap_fp16_exp(map_fp16_exp(v)) >= v
   static unsigned short map_fp16_exp(float v) {
     // float has exponents 10^{-44.85} .. 10^{38.53}
     int exp = (int)ceil(::log(v) / ::log(FP16_BASE)) + SHRT_UMAX / 2;
     if (exp < 0 || exp > SHRT_UMAX) {
       fprintf(stderr,"Error in map_fp16_exp(%g,%d)\n",v,exp);
       exit(3);
     }

     return (unsigned short)exp;
   }
  
   template<typename OPT>
     static void read_floats_fp16(char* & ptr, OPT* out, int64_t n, int nsc) {

     int64_t nsites = n / nsc;
     if (n % nsc) {
       fprintf(stderr,"Invalid size in write_floats_fp16\n");
       exit(4);
     }

     unsigned short* in = (unsigned short*)ptr;
     ptr += 2*(n+nsites);

     // do for each site
     for (int64_t site = 0;site<nsites;site++) {

       OPT* ev = &out[site*nsc];

       unsigned short* bptr = &in[site*(nsc + 1)];

       unsigned short exp = *bptr++;
       OPT max = unmap_fp16_exp(exp);
       OPT min = -max;

       for (int i=0;i<nsc;i++) {
	 ev[i] = fp_unmap( *bptr++, min, max, SHRT_UMAX );
       }

     }

   }

   template<typename OPT>
     static void write_floats_fp16(FILE* f, std::vector<char>& fbuf, uint32_t& crc, OPT* in, int64_t n, int nsc) {

     int64_t nsites = n / nsc;
     if (n % nsc) {
       fprintf(stderr,"Invalid size in write_floats_fp16\n");
       exit(4);
     }

     unsigned short* buf = (unsigned short*)malloc( sizeof(short) * (n + nsites) );
     if (!buf) {
       fprintf(stderr,"Out of mem\n");
       exit(1);
     }

     // do for each site
#pragma omp parallel for
     for (int64_t site = 0;site<nsites;site++) {

       OPT* ev = &in[site*nsc];

       unsigned short* bptr = &buf[site*(nsc + 1)];

       OPT max = fabs(ev[0]);
       OPT min;

       for (int i=0;i<nsc;i++) {
	 if (fabs(ev[i]) > max)
	   max = fabs(ev[i]);
       }

       unsigned short exp = map_fp16_exp(max);
       max = unmap_fp16_exp(exp);
       min = -max;

       *bptr++ = exp;

       for (int i=0;i<nsc;i++) {
	 int val = fp_map( ev[i], min, max, SHRT_UMAX );
	 if (val < 0 || val > SHRT_UMAX) {
	   fprintf(stderr,"Assert failed: val = %d (%d), ev[i] = %.15g, max = %.15g, exp = %d\n",val,SHRT_UMAX,ev[i],max,(int)exp);
	   exit(48);
	 }
	 *bptr++ = (unsigned short)val;
       }

     }

     write_bytes(buf,sizeof(short)*(n + nsites),f,fbuf,crc);

     free(buf);
   }

   template<typename Field,typename CoarseField>
     static bool read_compressed_vectors(const char* dir,BlockProjector<Field>& pr,BasisFieldVector<CoarseField>& coef, int ngroups = 1) {

     const BasisFieldVector<Field>& basis = pr._evec;
     GridBase* _grid = basis._v[0]._grid;

     // for error messages
     char hostname[1024];
     gethostname(hostname, 1024);

     std::cout << GridLogMessage << "Ready on host " << hostname << " with " << ngroups << " reader groups" << std::endl;

     // first read metadata
     char buf[4096];
     sprintf(buf,"%s/metadata.txt",dir);

     std::vector<uint32_t> s,b,nb,nn,crc32;
     s.resize(5); b.resize(5); nb.resize(5); nn.resize(5);
     uint32_t neig, nkeep, nkeep_single, blocks, _FP16_COEF_EXP_SHARE_FLOATS;
     uint32_t nprocessors = 1;

     FILE* f = 0;
     uint32_t status = 0;
     if (_grid->IsBoss()) {
       f = fopen(buf,"rb");
       status=f ? 1 : 0;
     }
     _grid->GlobalSum(status);
     std::cout << GridLogMessage << "Read params status " << status << std::endl;

     if (!status) {
       return false;
     }

#define _IRL_READ_INT(buf,p) if (f) { assert(fscanf(f,buf,p)==1); } else { *(p) = 0; } _grid->GlobalSum(*(p));

     for (int i=0;i<5;i++) {
       sprintf(buf,"s[%d] = %%d\n",i);
       _IRL_READ_INT(buf,&s[(i+1)%5]);
     }
     for (int i=0;i<5;i++) {
       sprintf(buf,"b[%d] = %%d\n",i);
       _IRL_READ_INT(buf,&b[(i+1)%5]);
     }
     for (int i=0;i<5;i++) {
       sprintf(buf,"nb[%d] = %%d\n",i);
       _IRL_READ_INT(buf,&nb[(i+1)%5]);
     }
     _IRL_READ_INT("neig = %d\n",&neig);
     _IRL_READ_INT("nkeep = %d\n",&nkeep);
     _IRL_READ_INT("nkeep_single = %d\n",&nkeep_single);
     _IRL_READ_INT("blocks = %d\n",&blocks);
     _IRL_READ_INT("FP16_COEF_EXP_SHARE_FLOATS = %d\n",&_FP16_COEF_EXP_SHARE_FLOATS);

     for (int i=0;i<5;i++) {
       assert(_grid->FullDimensions()[i] % s[i] == 0);
       nn[i] = _grid->FullDimensions()[i] / s[i];
       nprocessors *= nn[i];
     }

     std::cout << GridLogMessage << "Reading data that was generated on node-layout " << nn << std::endl;

     crc32.resize(nprocessors);
     for (int i =0;i<nprocessors;i++) {
       sprintf(buf,"crc32[%d] = %%X\n",i);
       _IRL_READ_INT(buf,&crc32[i]);
     }

#undef _IRL_READ_INT

     if (f)
       fclose(f);

     // allocate memory
     assert(std::equal(pr._bgrid._bs.begin(),pr._bgrid._bs.end(),b.begin())); // needs to be the same for now
     assert(pr._evec.size() == nkeep); // needs to be the same since this is compile-time fixed
     coef.resize(neig);


     // now get read geometry
     std::map<int, std::vector<int> > slots;
     std::vector<int> slot_lvol, lvol;
     int64_t slot_lsites;
     int ntotal;
     std::vector<int> _nn(nn.begin(),nn.end());
     get_read_geometry(_grid,_nn,
		       slots,slot_lvol,lvol,slot_lsites,
		       ntotal);
     int _nd = (int)lvol.size();

     // types
     typedef typename Field::scalar_type Coeff_t;
     typedef typename CoarseField::scalar_type CoeffCoarse_t;

     // slot layout
     int nperdir = ntotal / 32;
     if (nperdir < 1)
       nperdir=1;

     // add read groups
     for (int ngroup=0;ngroup<ngroups;ngroup++) {

       bool action = _grid->ThisRank() % ngroups == ngroup;
       
       std::cout << GridLogMessage << "Reading in group " << ngroup << " / " << ngroups << std::endl;

       // load all necessary slots and store them appropriately
       for (auto sl=slots.begin();sl!=slots.end();sl++) {

	 std::vector<int>& idx = sl->second;
	 int slot = sl->first;
	 std::vector<float> rdata;
       
	 char buf[4096];
       
	 if (action) {
	   // load one slot vector
	   sprintf(buf,"%s/%2.2d/%10.10d.compressed",dir,slot/nperdir,slot);
	   f = fopen(buf,"rb");
	   if (!f) {
	     fprintf(stderr,"Node %s cannot read %s\n",hostname,buf); fflush(stderr);
	     return false;
	   }
	 }

	 uint32_t crc = 0x0;
	 off_t size;
	   
	 GridStopWatch gsw;
	 _grid->Barrier();
	 gsw.Start();
	 
	 std::vector<char> raw_in(0);
	 if (action) {
	   fseeko(f,0,SEEK_END);
	   size = ftello(f);
	   fseeko(f,0,SEEK_SET);

	   raw_in.resize(size);
	   assert(fread(&raw_in[0],size,1,f) == 1);
	 }

	 _grid->Barrier();
	 gsw.Stop();

	 RealD totalGB = (RealD)size / 1024./1024./1024 * _grid->_Nprocessors;
	 RealD seconds = gsw.useconds() / 1e6;

	 if (action) {
	   std::cout << GridLogMessage << "[" << slot << "]  Read " << totalGB << " GB of compressed data at " << totalGB/seconds << " GB/s" << std::endl;
	   
	   uint32_t crc_comp = crc32_threaded((unsigned char*)&raw_in[0],size,0);
	   
	   if (crc_comp != crc32[slot]) {
	     std::cout << "Node " << hostname << " found crc mismatch for file " << buf << " (" << std::hex << crc_comp << " vs " << crc32[slot] << std::dec << ")" << std::endl;
	     std::cout << "Byte size: " << size << std::endl;
	   }

	   assert(crc_comp == crc32[slot]);
	 }

	 _grid->Barrier();
       
	 if (action) {
	   fclose(f);
	 }

	 char* ptr = &raw_in[0];
	 
	 GridStopWatch gsw2;
	 gsw2.Start();
	 if (action) {
	   int nsingleCap = nkeep_single;
	   if (pr._evec.size() < nsingleCap)
	     nsingleCap = pr._evec.size();
	   
	   int _cf_block_size = slot_lsites * 12 / 2 / blocks;
	   
#define FP_16_SIZE(a,b)  (( (a) + (a/b) )*2)
	   
	   // first read single precision basis vectors
#pragma omp parallel
	   {
	     std::vector<float> buf(_cf_block_size * 2);
#pragma omp for
	     for (int nb=0;nb<blocks;nb++) {
	       for (int i=0;i<nsingleCap;i++) {
		 char* lptr = ptr + buf.size()*(i + nsingleCap*nb)*4;
		 read_floats(lptr, &buf[0], buf.size() );
		 int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
		 if (mnb != -1)
		   pr._bgrid.pokeBlockOfVectorCanonical(mnb,pr._evec._v[i],buf);
	       }
	     }
	     
#pragma omp barrier
#pragma omp single
	     {
	       ptr = ptr + buf.size()*nsingleCap*blocks*4;
	     }
	     
	   }
	   
	   // TODO: at this point I should add a checksum test for block_sp(nb,v,v) for all blocks, then I would know that the mapping
	   // to blocks is OK at this point; after that ...
	   
	   // then read fixed precision basis vectors
#pragma omp parallel
	   {
	     std::vector<float> buf(_cf_block_size * 2);
#pragma omp for
	     for (int nb=0;nb<blocks;nb++) {
	       for (int i=nsingleCap;i<(int)pr._evec.size();i++) {
		 char* lptr = ptr + FP_16_SIZE( buf.size(), 24 )*((i-nsingleCap) + (pr._evec.size() - nsingleCap)*nb);
		 read_floats_fp16(lptr, &buf[0], buf.size(), 24);
		 int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
		 if (mnb != -1)
		   pr._bgrid.pokeBlockOfVectorCanonical(mnb,pr._evec._v[i],buf);
	       }
	     }
	     
#pragma omp barrier
#pragma omp single
	     {
	       ptr = ptr + FP_16_SIZE( buf.size()*(pr._evec.size() - nsingleCap)*blocks, 24 );
	     }
	   }
	   
#pragma omp parallel
	   {
	     std::vector<float> buf1(nkeep_single*2);
	     std::vector<float> buf2((nkeep - nkeep_single)*2);
	     
#pragma omp for
	     for (int j=0;j<(int)coef.size();j++)
	       for (int nb=0;nb<blocks;nb++) {
		 // get local coordinate on coarse grid
		 int ii,oi;
		 int mnb = pr._bgrid.globalToLocalCanonicalBlock(slot,_nn,nb);
		 if (mnb != -1)
		   canonical_block_to_coarse_coordinates(coef._v[0]._grid,mnb,ii,oi);
		 
		 char* lptr = ptr + (4*buf1.size() + FP_16_SIZE(buf2.size(), _FP16_COEF_EXP_SHARE_FLOATS))*(nb + j*blocks);
		 int l;
		 read_floats(lptr, &buf1[0], buf1.size() );
		 if (mnb != -1) {
		   for (l=0;l<nkeep_single;l++) {
		     ((CoeffCoarse_t*)&coef._v[j]._odata[oi]._internal._internal[l])[ii] = CoeffCoarse_t(buf1[2*l+0],buf1[2*l+1]);
		   }
		 }
		 read_floats_fp16(lptr, &buf2[0], buf2.size(), _FP16_COEF_EXP_SHARE_FLOATS);
		 if (mnb != -1) {
		   for (l=nkeep_single;l<nkeep;l++) {
		     ((CoeffCoarse_t*)&coef._v[j]._odata[oi]._internal._internal[l])[ii] = CoeffCoarse_t(buf2[2*(l-nkeep_single)+0],buf2[2*(l-nkeep_single)+1]);
		   }
		 }
		 
	       }
	   }
	   
	   // set checkerboard
	   for (int i=0;i<(int)pr._evec.size();i++)
	     pr._evec._v[i].checkerboard = Odd;
	 
	   gsw2.Stop();
	   seconds=gsw2.useconds()/1e6;
	   std::cout << GridLogMessage << "Processed " << totalGB << " GB of compressed data at " << totalGB/seconds << " GB/s" << std::endl;
	 }
       }
     }
#undef FP_16_SIZE
     return true;

   }
 
   static bool DirectoryExists(const char *path) {
     struct stat info;
     return ((stat( path, &info ) == 0) && (info.st_mode & S_IFDIR));
   }
 
   static void conditionalMkDir(const char* path) {
     if (!DirectoryExists(path))
       mkdir(path,ACCESSPERMS);
   }
 
   template<typename Field,typename CoarseField>
   static void write_compressed_vectors(const char* dir,const BlockProjector<Field>& pr,
					const BasisFieldVector<CoarseField>& coef,
					int nsingle,int writer_nodes = 0) {

     GridStopWatch gsw;
     
     const BasisFieldVector<Field>& basis = pr._evec;
     GridBase* _grid = basis._v[0]._grid;
     std::vector<int> _l = _grid->FullDimensions();
     for (int i=0;i<(int)_l.size();i++)
       _l[i] /= _grid->_processors[i];

     _grid->Barrier();
     gsw.Start();

     char buf[4096];
     
     // Making the directories is somewhat tricky.
     // If we run on a joint filesystem we would just 
     // have the boss create the directories and then
     // have a barrier.  We also want to be able to run
     // on local /scratch, so potentially all nodes need
     // to create their own directories.  So do the following
     // for now.
     for (int j=0;j<_grid->_Nprocessors;j++) {
       if (j == _grid->ThisRank()) {
	 conditionalMkDir(dir);
	 for (int i=0;i<32;i++) {
	   sprintf(buf,"%s/%2.2d",dir,i);
	   conditionalMkDir(buf);
	 }       
	 _grid->Barrier(); // make sure directories are ready
       }
     }
     

     typedef typename Field::scalar_type Coeff_t;
     typedef typename CoarseField::scalar_type CoeffCoarse_t;
     
     int nperdir = _grid->_Nprocessors / 32;
     if (nperdir < 1)
       nperdir=1;

     int slot;
     Lexicographic::IndexFromCoor(_grid->_processor_coor,slot,_grid->_processors);

     int64_t off = 0x0;
     uint32_t crc = 0x0;
     if (writer_nodes < 1)
       writer_nodes = _grid->_Nprocessors;
     int groups = _grid->_Nprocessors / writer_nodes;
     if (groups<1)
       groups = 1;

     std::cout << GridLogMessage << " Write " << dir << " nodes = " << writer_nodes << std::endl;

     for (int group=0;group<groups;group++) {
       _grid->Barrier();
       if (_grid->ThisRank() % groups == group) {
	 
	 sprintf(buf,"%s/%2.2d/%10.10d.compressed",dir,slot/nperdir,slot);
	 FILE* f = fopen(buf,"wb");
	 assert(f);

	 //buffer does not seem to help
	 //assert(!setvbuf ( f , NULL , _IOFBF , 1024*1024*2 ));
	 	 
	 int nsingleCap = nsingle;
	 if (pr._evec.size() < nsingleCap)
	   nsingleCap = pr._evec.size();
	 
	 GridStopWatch gsw1,gsw2,gsw3,gsw4,gsw5;

	 gsw1.Start();

	 std::vector<char> fbuf;
	 fbuf.reserve( 1024 * 1024 * 8 );

	 // first write single precision basis vectors
	 for (int nb=0;nb<pr._bgrid._blocks;nb++) {
	   for (int i=0;i<nsingleCap;i++) {
	     std::vector<float> buf;
	     pr._bgrid.peekBlockOfVectorCanonical(nb,pr._evec._v[i],buf);
	     
#if 0
	     {
	       RealD nrm = 0.0;
	       for (int j=0;j<(int)buf.size();j++)
		 nrm += buf[j]*buf[j];
	       std::cout << GridLogMessage << "Norm: " << nrm << std::endl;
	     }
#endif
	     write_floats(f,fbuf,crc, &buf[0], buf.size() );
	   }
	 }

	 gsw1.Stop();
	 gsw2.Start();

	 // then write fixed precision basis vectors
	 for (int nb=0;nb<pr._bgrid._blocks;nb++) {
	   for (int i=nsingleCap;i<(int)pr._evec.size();i++) {
	     std::vector<float> buf;
	     pr._bgrid.peekBlockOfVectorCanonical(nb,pr._evec._v[i],buf);
	     write_floats_fp16(f,fbuf,crc, &buf[0], buf.size(), 24);
	   }
	 }

	 gsw2.Stop();
	 assert(coef._v[0]._grid->_isites*coef._v[0]._grid->_osites == pr._bgrid._blocks);

	 gsw3.Start();
	 for (int j=0;j<(int)coef.size();j++) {

	   int64_t size1 = nsingleCap*2;
	   int64_t size2 = 2*(pr._evec.size()-nsingleCap);
	   int64_t size = size1;
	   if (size2>size)
	     size=size2;
	   std::vector<float> buf(size);

	   //RealD nrmTest = 0.0;
	   for (int nb=0;nb<pr._bgrid._blocks;nb++) {
	     // get local coordinate on coarse grid
	     int ii, oi;
	     canonical_block_to_coarse_coordinates(coef._v[0]._grid,nb,ii,oi);
	     
	     gsw4.Start();
	     gsw5.Start();
	     for (int l=0;l<nsingleCap;l++) {
	       auto res = ((CoeffCoarse_t*)&coef._v[j]._odata[oi]._internal._internal[l])[ii];
	       buf[2*l+0] = res.real();
	       buf[2*l+1] = res.imag();
	       //nrmTest += res.real() * res.real() + res.imag() * res.imag();
	     }
	     gsw5.Stop();
	     write_floats(f,fbuf,crc, &buf[0], size1 );
	     gsw4.Stop();
	     for (int l=nsingleCap;l<(int)pr._evec.size();l++) {
	       auto res = ((CoeffCoarse_t*)&coef._v[j]._odata[oi]._internal._internal[l])[ii];
	       buf[2*(l-nsingleCap)+0] = res.real();
	       buf[2*(l-nsingleCap)+1] = res.imag();
	       //nrmTest += res.real() * res.real() + res.imag() * res.imag();
	     }
	     write_floats_fp16(f,fbuf,crc, &buf[0], size2, FP16_COEF_EXP_SHARE_FLOATS);
	   }

	   //_grid->GlobalSum(nrmTest);
	   //std::cout << GridLogMessage << "Test norm: " << nrmTest << std::endl;
	 }
	 gsw3.Stop();

	 flush_bytes(f,fbuf);

	 off = ftello(f);	 
	 fclose(f);

	 std::cout<<GridLogMessage << "Timing: write single basis in " << gsw1.Elapsed() << " and fp16 in " << gsw2.Elapsed() << " and coefficients in " << gsw3.Elapsed() << " (" <<
	   gsw4.Elapsed() << ", " << gsw5.Elapsed() << ")" << std::endl;
       }
     }
	 
     _grid->Barrier();
     gsw.Stop();
     
     RealD totalGB = (RealD)off / 1024./1024./1024 * _grid->_Nprocessors;
     RealD seconds = gsw.useconds() / 1e6;
     std::cout << GridLogMessage << "Write " << totalGB << " GB of compressed data at " << totalGB/seconds << " GB/s in " << seconds << " s" << std::endl;

     // gather crcs
     std::vector<uint32_t> crcs(_grid->_Nprocessors);
     for (int i=0;i<_grid->_Nprocessors;i++) {
       crcs[i] = 0x0;
     }
     crcs[slot] = crc;
     for (int i=0;i<_grid->_Nprocessors;i++) {
       _grid->GlobalSum(crcs[i]);
     }
     
     if (_grid->IsBoss()) {
       sprintf(buf,"%s/metadata.txt",dir);
       FILE* f = fopen(buf,"wb");
       assert(f);
       for (int i=0;i<5;i++)
	 fprintf(f,"s[%d] = %d\n",i,_grid->FullDimensions()[(i+1)%5] / _grid->_processors[(i+1)%5]);
       for (int i=0;i<5;i++)
	 fprintf(f,"b[%d] = %d\n",i,pr._bgrid._bs[(i+1)%5]);
       for (int i=0;i<5;i++)
	 fprintf(f,"nb[%d] = %d\n",i,pr._bgrid._nb[(i+1)%5]);
       fprintf(f,"neig = %d\n",(int)coef.size());
       fprintf(f,"nkeep = %d\n",(int)pr._evec.size());
       fprintf(f,"nkeep_single = %d\n",nsingle);
       fprintf(f,"blocks = %d\n",pr._bgrid._blocks);
       fprintf(f,"FP16_COEF_EXP_SHARE_FLOATS = %d\n",FP16_COEF_EXP_SHARE_FLOATS);
       for (int i =0;i<_grid->_Nprocessors;i++)
	 fprintf(f,"crc32[%d] = %X\n",i,crcs[i]);
       fclose(f);
     }

   }

   template<typename Field>
     static void write_argonne(const BasisFieldVector<Field>& ret,const char* dir) {
     
     GridBase* _grid = ret._v[0]._grid;
     std::vector<int> _l = _grid->FullDimensions();
     for (int i=0;i<(int)_l.size();i++)
       _l[i] /= _grid->_processors[i];

     char buf[4096];
     
     if (_grid->IsBoss()) {
       mkdir(dir,ACCESSPERMS);
       
       for (int i=0;i<32;i++) {
	 sprintf(buf,"%s/%2.2d",dir,i);
	 mkdir(buf,ACCESSPERMS);
       }
     }
     
     _grid->Barrier(); // make sure directories are ready

     
     int nperdir = _grid->_Nprocessors / 32;
     if (nperdir < 1)
       nperdir=1;
     std::cout << GridLogMessage << " Write " << dir << " nodes = " << _grid->_Nprocessors << std::endl;
     
     int slot;
     Lexicographic::IndexFromCoor(_grid->_processor_coor,slot,_grid->_processors);
     //printf("Slot: %d <> %d\n",slot, _grid->ThisRank());
     
     sprintf(buf,"%s/%2.2d/%10.10d",dir,slot/nperdir,slot);
     FILE* f = fopen(buf,"wb");
     assert(f);
     
     int N = (int)ret._v.size();
     uint32_t crc = 0x0;
     int64_t cf_size = _grid->oSites()*_grid->iSites()*12;
     std::vector< float > rdata(cf_size*2);
     
     GridStopWatch gsw1,gsw2;
     
     for (int i=0;i<N;i++) {

       if (ret._v[i].checkerboard != Odd)
	 continue;

       // create buffer and put data in argonne format in there
       std::vector<int> coor(_l.size());
       for (coor[1] = 0;coor[1]<_l[1];coor[1]++) {
	 for (coor[2] = 0;coor[2]<_l[2];coor[2]++) {
	   for (coor[3] = 0;coor[3]<_l[3];coor[3]++) {
	     for (coor[4] = 0;coor[4]<_l[4];coor[4]++) {
	       for (coor[0] = 0;coor[0]<_l[0];coor[0]++) {
		 
		 if ((coor[1]+coor[2]+coor[3]+coor[4]) % 2 == 1) {
		   // peek
		   iScalar<iVector<iVector<ComplexF, 3>, 4> > sc;
		   peekLocalSite(sc,ret._v[i],coor);
		   for (int s=0;s<4;s++)
		     for (int c=0;c<3;c++)
		       *(std::complex<float>*)&rdata[get_bfm_index(&coor[0],c+s*3, &_l[0] )] = sc()(s)(c);
		 }
	       }
	     }
	   }
	 }
       }
       
       // endian flip
       for (int i=0;i<cf_size*2;i++) {
	 char* c = (char*)&rdata[i];
	 char tmp; int j;
	 for (j=0;j<2;j++) {
	   tmp = c[j]; c[j] = c[3-j]; c[3-j] = tmp;
	 }
       }
       
       // create crc of buffer
       gsw1.Start();
       crc = crc32_threaded((unsigned char*)&rdata[0],cf_size*2*4,crc);    
       gsw1.Stop();
       
       // write out
       gsw2.Start();
       assert(fwrite(&rdata[0],cf_size*2*4,1,f)==1);
       gsw2.Stop();
       
     }
     
     fclose(f);
     
     
     // gather crc's and write out
     std::vector<uint32_t> crcs(_grid->_Nprocessors);
     for (int i=0;i<_grid->_Nprocessors;i++) {
       crcs[i] = 0x0;
     }
     crcs[slot] = crc;
     for (int i=0;i<_grid->_Nprocessors;i++) {
       _grid->GlobalSum(crcs[i]);
     }
     
     if (_grid->IsBoss()) {
       sprintf(buf,"%s/checksums.txt",dir);
       FILE* f = fopen(buf,"wt");
       assert(f);
       fprintf(f,"00000000\n\n");
       for (int i =0;i<_grid->_Nprocessors;i++)
	 fprintf(f,"%X\n",crcs[i]);
       fclose(f);
       
       sprintf(buf,"%s/nodes.txt",dir);
       f = fopen(buf,"wt");
       assert(f);
       for (int i =0;i<(int)_grid->_processors.size();i++)
	 fprintf(f,"%d\n",_grid->_processors[i]);
       fclose(f);
     }
     
     
     std::cout << GridLogMessage << "Writing slot " << slot << " with "
	       << N << " vectors in "
	       << gsw2.Elapsed() << " at " 
	       << ( (double)cf_size*2*4 * N / 1024./1024./1024. / gsw2.useconds()*1000.*1000. )
	       << " GB/s  with crc computed at "
	       << ( (double)cf_size*2*4 * N / 1024./1024./1024. / gsw1.useconds()*1000.*1000. )
	       << " GB/s "
	       << std::endl;
     
     _grid->Barrier();
     std::cout << GridLogMessage  << "Writing complete" << std::endl;
     
   }
 }

}
