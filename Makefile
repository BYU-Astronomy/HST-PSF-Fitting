INCLUDES = -I${CONDA_PREFIX}/include
CXXFLAGS = -O3 -g
XSIMDFLAGS = -mavx2 -ffast-math -DXTENSOR_USE_XSIMD

.PHONY: clean conda-setup-build conda-setup-run

all: PSF_prep Single_Fitting Binary_Fitting Second_Binary_Fitting Secondary_Binary_Fitting

PSF_prep: PSF_prep.cpp
	g++ ${CXXFLAGS} ${XSIMDFLAGS} ${INCLUDES} -o $@ $^

Single_Fitting: Single_Fitting.cpp
	g++ ${CXXFLAGS} ${XSIMDFLAGS} ${INCLUDES} -o $@ $^

Binary_Fitting: Binary_Fitting.cpp
	g++ ${CXXFLAGS} ${XSIMDFLAGS} ${INCLUDES} -o $@ $^

Second_Binary_Fitting: Second_Binary_Fitting.cpp
	g++ ${CXXFLAGS} ${XSIMDFLAGS} ${INCLUDES} -o $@ $^

Secondary_Binary_Fitting: Secondary_Binary_Fitting.cpp
	g++ ${CXXFLAGS} ${XSIMDFLAGS} ${INCLUDES} -o $@ $^

clean:
	rm -f PSF_prep
	rm -f Single_Fitting
	rm -f Binary_Fitting
	rm -f Second_Binary_Fitting
	rm -f Secondary_Binry_Fitting


conda-setup-build:
	conda install -c conda-forge xtensor xsimd numpy astropy photutils

conda-setup-run:
	conda install -c conda-forge numpy astropy photutils
