CC=g++
CFLAGS=`root-config --evelibs --cflags`

showGeometry:
	$(CC) -o showGeometry showGeometry.cpp $(CFLAGS)

.PHONY: clean

clean:
	rm -f showGeometry
