.PHONY: all
all: visibility corrsift

.PHONY: visibility
visibility:
	$(MAKE) -C rrn/visibility visibilitylibomp

.PHONY: corrsift
corrsift:
	$(MAKE) -C rrn/corrsift libcorrsift.so

.PHONY: cleanall
cleanall:
	$(MAKE) -C rrn/visibility clean
	$(MAKE) -C rrn/corrsift clean
