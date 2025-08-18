format:
	ruff check --exit-zero --fix .
	ruff format .

run:
	rm -rf main_output
	python main.py

view:
	pysph view main_output
