format:
	ruff check --exit-zero --fix .
	ruff format .

run:
	rm -rf /home/isadorarodrigues/sph/main_output
	python main.py

view:
	pysph view /home/isadorarodrigues/sph/main_output
