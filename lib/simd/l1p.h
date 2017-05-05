#pragma once
namespace Grid {
// L1p optimisation 
inline void bgq_l1p_optimisation(int mode)
{
#ifdef QPX
#undef L1P_CFG_PF_USR
#define L1P_CFG_PF_USR  (0x3fde8000108ll)   /*  (64 bit reg, 23 bits wide, user/unpriv) */

  uint64_t cfg_pf_usr;
  if ( mode ) { 
    cfg_pf_usr =
        L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_adaptive_throttle(0xF) ;
    //    if ( sizeof(Float) == sizeof(double) ) {
      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(2)| L1P_CFG_PF_USR_dfetch_max_footprint(3)   ;
      //    } else {
      //      cfg_pf_usr |=  L1P_CFG_PF_USR_dfetch_depth(1)| L1P_CFG_PF_USR_dfetch_max_footprint(2)   ;
      //    }
  } else { 
    cfg_pf_usr = L1P_CFG_PF_USR_dfetch_depth(1)
      | L1P_CFG_PF_USR_dfetch_max_footprint(2)   
      | L1P_CFG_PF_USR_ifetch_depth(0)       
      | L1P_CFG_PF_USR_ifetch_max_footprint(1)   
      | L1P_CFG_PF_USR_pf_stream_est_on_dcbt 
      | L1P_CFG_PF_USR_pf_stream_establish_enable
      | L1P_CFG_PF_USR_pf_stream_optimistic
      | L1P_CFG_PF_USR_pf_stream_prefetch_enable;
  }
  *((uint64_t *)L1P_CFG_PF_USR) = cfg_pf_usr;
#endif
}
}
