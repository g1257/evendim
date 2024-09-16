# This file is to be included by ../Config.make
#XACC_INSTALL_PATH=$(HOME)/.xacc

CPPFLAGS += -I $(XACC_INSTALL_PATH)/include  -I $(XACC_INSTALL_PATH)/include/xacc -I$(XACC_INSTALL_PATH)/include/cppmicroservices3 -I$(XACC_INSTALL_PATH)/include/quantum/gate -DUSE_XACC
CPPFLAGS += -I $(XACC_INSTALL_PATH)/../xacc/xacc/accelerator/remote

LDFLAGS  += -L$(XACC_INSTALL_PATH)/lib -lxacc-pauli -lCppMicroServices -lxacc

