FC         = nvfortran
FFLAGS     = -O2 -module $(BUILDIR)
CUDA_FLAGS = -cuda

BUILDIR = bin
TARGET  = vec_mat_add

OBJS = utils_module.o vec_add_cpu_module.o mat_add_cpu_module.o vec_add_module.o mat_add_module.o main.o

all: $(BUILDIR)/$(TARGET)

$(BUILDIR):
	mkdir -p $(BUILDIR)

$(BUILDIR)/$(TARGET): $(addprefix $(BUILDIR)/,$(OBJS))
	$(FC) $(FFLAGS) $(CUDA_FLAGS) -o $@ $^

$(BUILDIR)/utils_module.o: utils_module.f90 | $(BUILDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDIR)/vec_add_cpu_module.o: vec_add_cpu.f90 | $(BUILDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDIR)/mat_add_cpu_module.o: mat_add_cpu.f90 | $(BUILDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDIR)/vec_add_module.o: vec_add_module.cuf | $(BUILDIR)
	$(FC) $(FFLAGS) $(CUDA_FLAGS) -c $< -o $@

$(BUILDIR)/mat_add_module.o: mat_add_module.cuf | $(BUILDIR)
	$(FC) $(FFLAGS) $(CUDA_FLAGS) -c $< -o $@

$(BUILDIR)/main.o: main.f90 $(BUILDIR)/utils_module.o $(BUILDIR)/vec_add_cpu_module.o $(BUILDIR)/mat_add_cpu_module.o $(BUILDIR)/vec_add_module.o $(BUILDIR)/mat_add_module.o
	$(FC) $(FFLAGS) $(CUDA_FLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDIR)/,$(OBJS)) $(BUILDIR)/*.mod $(BUILDIR)/$(TARGET)
