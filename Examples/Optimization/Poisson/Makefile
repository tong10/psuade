target:
	cc -O -o smmm_fine_sim smmm_fine_sim.c -lm
	cc -O -o smmm_coarse_sim smmm_coarse_sim.c -lm
	cc -O -o cobyla_simulator cobyla_simulator.c -lm

clean:
	@rm -f *.o auxData
	@rm -f cobyla_simulator smmm_coarse_sim smmm_fine_sim 
	@rm -f psuadeData* psuadeApps* psuadeOpt*

