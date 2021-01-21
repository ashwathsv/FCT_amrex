AMREX_HOME ?= /Users/ashwathsv/amrex
PRECISION  = DOUBLE
# TINY_PROFILE = TRUE
# PROFILE      = FALSE
# AMREX_ENABLE_TINY_PROFILE = ON
# AMREX_ENABLE_BASE_PROFILE = OFF

DEBUG      = TRUE
# DEBUG      = FALSE

# DIM        = 2
DIM       = 3

COMP	   = gnu
CXXSTD     = c++11
INCLUDE_LOCATIONS += /usr/local/opt/llvm/include
LIBRARY_LOCATIONS += /usr/local/opt/llvm/lib
LIBRARIES += -lomp

USE_MPI    = TRUE
USE_OMP    = FALSE
# USE_OMP    = TRUE

Bpack   := ./Make.package 
Blocs   := . 

include ../Make.Adv
