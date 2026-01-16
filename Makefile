format:
	ruff check --exit-zero --fix .
	ruff format .

run:
	rm -rf main_output
	python main.py

view:
	pysph view main_output

paraview:
	pysph dump_vtk main_output -d /Users/isa/Documents/sph_results

run_view:
	make run
	make view
