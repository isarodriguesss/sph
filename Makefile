format:
	ruff check --exit-zero --fix .
	ruff format .

run:
	rm -rf main_output
	python main.py
	python plot.py

view:
	pysph view main_output

paraview:
	rm -rf /Users/isa/Developer/sph_results/paraview
	pysph dump_vtk main_output -d /Users/isa/Developer/sph_results/paraview

run_view:
	make run
	make view

plot:
	python plot.py
