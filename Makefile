
all : overall prior_pop example coverage_sim

dirs :
	mkdir -p figures
	mkdir -p saved_results
	mkdir -p tables

overall : dirs
	R < overall.R > saved_results/overall.out.txt --no-save
	R < overall_suest.R > saved_results/overall_suest.out.txt --no-save

prior_pop : dirs
	R < by_popularity.R > saved_results/by_popularity.out.txt --no-save

example : dirs
	R < example_domain_plots.R > saved_results/example_domain_plots.out.txt --no-save

coverage_sim : dirs
	R < coverage_sim.R > saved_results/coverage_sim.out.txt --no-save

clean :
	rm -fR saved_results
	rm -fR figures
	rm -fR tables
	rm -f Rplots.pdf
