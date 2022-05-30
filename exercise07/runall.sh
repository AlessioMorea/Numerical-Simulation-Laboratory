cd MDNVE_MCNVT
bash clean.sh
make
./NVE_NVT.exe gas
./NVE_NVT.exe liquid
./NVE_NVT.exe solid
./NVE_NVT.exe gasNotEq
./NVE_NVT.exe liquidNotEq
./NVE_NVT.exe solidNotEq
./NVE_NVT.exe gasNVE
./NVE_NVT.exe liquidNVE
./NVE_NVT.exe solidNVE