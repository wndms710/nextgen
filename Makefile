mpi ?= 1

include Makefile.in
include makefiles/Makefile.version

default:
		@echo "********************************************************************"
		@echo "*                                                                  *"
		@echo "* Copyright 1997-2020 Cascade Technologies Inc.                    *"
		@echo "* nextgen version" $(DOCS_VERSION)                                              
		@echo "*                                                                  *"
		@echo "********************************************************************"
		@echo "*                                                                  *"
		@echo "* The executables, source code and documentation (collectively the *"
		@echo "* 'Cascade Software') are licensed by Cascade Technologies for use *"
		@echo "* only under the terms of the licence agreement.                   *"
		@echo "*                                                                  *"
		@echo "* No part of the Cascade Software may be stored, reproduced, or    *"
		@echo "* transmitted in any form or by any means without the prior        *"
		@echo "* written permission of Cascade Technologies.                      *"
		@echo "*                                                                  *"
		@echo "********************************************************************"
		make -C ./src/surfer
		make -C ./src/stitch
		make -C ./src/charles
		make -C ./src/tools
		make -C ./src/ping
		make -C ./src/acoustics

tests:
		make -C ./src/stitch/tests
		make -C ./src/charles/tests 

clean:
		make -C ./src/core/common clean
		make -C ./src/core/encrypt clean
		make -C ./src/surfer clean
		make -C ./src/stitch clean
		make -C ./src/charles clean
		make -C ./src/tools clean
		make -C ./src/ping clean
		make -C ./src/chemistry clean
		make -C ./src/helmholtz_vof clean
		make -C ./src/vof clean
		make -C ./tests clean

#RELEASE_DIR?=./release
#release:
#	@echo "Release directory: $(RELEASE_DIR)"
#	@if [ -d "$(RELEASE_DIR)" ]; then echo "RELEASE_DIR exists, remove before re-running"; exit 1 ; fi
#	@echo "Building tools..."
#	make -C ./src/tools
#	@echo "Creating release directory structure..."
#	mkdir -p $(RELEASE_DIR)/bin
#	mkdir -p $(RELEASE_DIR)/src/core
#	mkdir $(RELEASE_DIR)/src/tools
#	mkdir $(RELEASE_DIR)/src/charles
#	mkdir $(RELEASE_DIR)/src/chemistry
#	@echo "Copying core libraries and headers..."
#	cp -r ./src/core/lib     $(RELEASE_DIR)/src/core
#	cp -r ./src/core/include $(RELEASE_DIR)/src/core
#	@echo "Copying solver libraries..."
#	cp -r ./src/tools/json_templates $(RELEASE_DIR)/src/tools
#	cp ./src/chemistry/*hpp  $(RELEASE_DIR)/src/chemistry
#	cp ./src/chemistry/libchemistry.a $(RELEASE_DIR)/src/chemistry
#	cp ./src/charles/*hpp $(RELEASE_DIR)/src/charles
#	cp ./src/charles/charles.cpp $(RELEASE_DIR)/src/charles
#	@echo "Copying release executables..."
#	cp -r ./bin $(RELEASE_DIR)
#	@echo "Copying makefiles..."
#	cp Makefile.in $(RELEASE_DIR)
#	cp makefiles/Makefile.release $(RELEASE_DIR)
#	@sed -i 's@CTI_RELEASE_HOME?=@CTI_RELEASE_HOME?=$(RELEASE_DIR)@' $(RELEASE_DIR)/Makefile.release
#	@echo "**************************REMINDER**************************"
#	@echo "createChemtable.exe and stitch.exe must be deployed manually"
#	@echo "************************************************************"


.PHONY: clean tests
