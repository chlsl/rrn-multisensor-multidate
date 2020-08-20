STRICT=# -fPIC -ansi -Wall -Wextra -Werror
IIO=-DNDEBUG -std=c99 -lpng -ltiff -ljpeg
OPT= -O3
OpenMP=

default:
	@echo "make clouds"
	@echo "or make cloudsomp"
	@echo "or make cloudslib.so"
	@echo "or make cloudslibomp"
	@echo "or <same commands but replace clouds by visibility>"

addopenmp:
ifeq ($(shell uname), Darwin)
OpenMP=-Xpreprocessor -fopenmp -lomp
else
OpenMP=-fopenmp
endif

clouds: clouds.c iio.c
	$(CC) -o $@ $^ $(OPT) $(OpenMP) $(STRICT) $(IIO) -lm

cloudslib.so: clouds.c iio.c
	$(CC) -shared -fPIC -o $@ $^ $(OPT) $(OpenMP) $(STRICT) $(IIO) -lm

cloudsomp: addopenmp
cloudsomp: clouds

# cloudslibomp: OpenMP=-Xpreprocessor -fopenmp -lomp # CHLSL: for mac os
cloudslibomp: addopenmp
cloudslibomp: cloudslib.so

visibility: visibility.c iio.c
	$(CC) $(OPT) $(OpenMP) $(STRICT) -o $@ $^ $(IIO) -lm

visibilitylib.so: visibility.c iio.c
	$(CC) -DASLIB -shared -fPIC -o $@ $^ $(OPT) $(OpenMP) $(STRICT) $(IIO) -lm

visibilitylibomp: addopenmp
visibilitylibomp: visibilitylib.so

visibilityomp: addopenmp
visibilityomp: visibility

clean:
	$(RM) clouds visibility *.o *.so

