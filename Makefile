format:
	ruff check --exit-zero --fix .
	ruff format .

run:
	rm -rf main_output
	python main.py
	python plots/plot.py

view:
	pysph view main_output

paraview:
	rm -rf /Users/isa/Developer/sph_results/paraview
	pysph dump_vtk main_output -d /Users/isa/Developer/sph_results/paraview

run_view:
	make run
	make view

plot:
	python plots/plot.py

RUN ?= runs/E11_t100
T ?= 50
fig:
	python plots/fig_tese.py $(RUN) --t $(T)
