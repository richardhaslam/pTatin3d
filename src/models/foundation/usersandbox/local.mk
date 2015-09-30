libptatin3dmodels-y.c += $(call thisdir, \
    userevaluator_empty.c \
    useroutput_materialpointstd_viewer.c \
   	)

TATIN_INC += -I$(abspath $(call thisdir,.))


