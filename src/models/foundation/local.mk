libptatin3dmodels-y.c += $(call thisdir, \
			foundation_reg.c \
			foundation_init.c \
			foundation_destroy.c \
            foundation_meshic.c \
            foundation_rheology.c \
            foundation_materialgeomic.c \
            foundation_mpeval.c \
            foundation_user.c \
            foundation_userfunctions_reg.c \
	)

TATIN_INC += -I$(abspath $(call thisdir,.))

include $(call incsubdirs,usersandbox)

