# Nothing
CPPFLAGS += -I$(HOME)/.xacc/include  -I $(HOME)/.xacc/include/xacc -I$(HOME)/.xacc/include/cppmicroservices4 -I$(HOME)/.xacc/include/quantum/gate -DUSE_XACC
LDFLAGS  += -L$(HOME)/.xacc/lib -lxacc

