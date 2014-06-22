all: c3reliability.py replicatedstar.py

c3reliability.py: c3reliability.sage prismmodel.py
	sage c3reliability.sage

replicatedstar.py: replicatedstar.sage prismmodel.py
	sage replicatedstar.sage

prismmodel.py: prismmodel.sage
	sage prismmodel.sage

.PHONY: clean

clean:
