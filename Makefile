.PHONY: all cmake-debug debug cmake-release release clean

BUILD_DIR	:= _build
RELEASE_DIR	:= $(BUILD_DIR)/release
DEBUG_DIR 	:= $(BUILD_DIR)/debug

CMAKE_FLAGS	:= -D CMAKE_C_COMPILER=$(CC)

all: release

cmake-debug:
	mkdir -p $(DEBUG_DIR)
	cd $(DEBUG_DIR) && cmake $(CMAKE_FLAGS) -D CMAKE_BUILD_TYPE=Debug ../..

debug: cmake-debug
	$(MAKE) -C $(DEBUG_DIR)

cmake-release:
	mkdir -p $(RELEASE_DIR)
	cd $(RELEASE_DIR) && cmake $(CMAKE_FLAGS) -D CMAKE_BUILD_TYPE=Release ../..

release: cmake-release
	$(MAKE) -C $(RELEASE_DIR)

clean:
	rm -rf $(BUILD_DIR)
