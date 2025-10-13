.PHONY: all clean

# 需要构建的目录
TARGET_DIRS = \
	utils/faSize \
	utils/faSplit \
	utils/faToTwoBit \
	hg/mouseStuff/axtChain \
	hg/mouseStuff/axtToMaf \
	hg/mouseStuff/chainNet \
	hg/mouseStuff/chainSort \
	hg/mouseStuff/netToAxt \
	hg/ratStuff/mafFilter \
	hg/ratStuff/mafSplit \
	hg/chainNet \
	hg/trfBig

# 默认构建全部
all:
	@for dir in $(TARGET_DIRS); do \
		echo "Building in $$dir..."; \
		( cd $$dir && $(MAKE) ); \
	done

# 清理全部
clean:
	@for dir in $(TARGET_DIRS); do \
		echo "Cleaning in $$dir..."; \
		( cd $$dir && $(MAKE) clean || true ); \
	done
