ifeq ($(strip $(analyzers/analyzers)),)
ALL_COMMONRULES += src_analyzers_analyzers_src
src_analyzers_analyzers_src_parent := analyzers/analyzers
src_analyzers_analyzers_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_analyzers_analyzers_src,src/analyzers/analyzers/src,LIBRARY))
analyzersanalyzers := self/analyzers/analyzers
analyzers/analyzers := analyzersanalyzers
analyzersanalyzers_files := $(patsubst src/analyzers/analyzers/src/%,%,$(wildcard $(foreach dir,src/analyzers/analyzers/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
analyzersanalyzers_BuildFile    := $(WORKINGDIR)/cache/bf/src/analyzers/analyzers/BuildFile
analyzersanalyzers_LOC_USE := self  FWCore/Framework DataFormats/Candidate FWCore/PluginManager FWCore/ParameterSet PhysicsTools/CandUtils PhysicsTools/UtilAlgos fastjet
analyzersanalyzers_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,analyzersanalyzers,analyzersanalyzers,$(SCRAMSTORENAME_LIB),src/analyzers/analyzers/src))
analyzersanalyzers_PACKAGE := self/src/analyzers/analyzers/src
ALL_PRODS += analyzersanalyzers
analyzersanalyzers_CLASS := LIBRARY
analyzers/analyzers_forbigobj+=analyzersanalyzers
analyzersanalyzers_INIT_FUNC        += $$(eval $$(call Library,analyzersanalyzers,src/analyzers/analyzers/src,src_analyzers_analyzers_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
endif
