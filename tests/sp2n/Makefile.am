.PHONY: check

include Make.inc

check: tests
	./Test_project_on_Sp
	./Test_sp2n_lie_gen
	./Test_Sp_start
	./Test_hmc_Sp_Wilson2ASFermionGauge --grid 4.4.4.4 --Trajectories 1 --Thermalizations 0 --StartingType ColdStart --assert=true
