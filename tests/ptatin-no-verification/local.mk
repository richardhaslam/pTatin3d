ptatin-tests-y.c += $(call thisdir, \
			test_cjson.c \
	)

TATIN_INC += -I$(abspath $(call thisdir,.))

