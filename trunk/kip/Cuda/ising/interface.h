struct Ising;

extern "C" void ising_cuda_init_device();
extern "C" Ising *ising_cuda_create(int len, int dim, float h, float T);
extern "C" void ising_cuda_destroy(Ising *ising);
extern "C" double ising_cuda_magnetization(Ising *ising);
extern "C" void ising_cuda_get_spins(Ising *ising, int *spins);
extern "C" void ising_cuda_set_spins(Ising *ising, int *spins);
extern "C" void ising_cuda_step(Ising *ising, int parity);
extern "C" void ising_cuda_quench(Ising *ising, float h, float T);
