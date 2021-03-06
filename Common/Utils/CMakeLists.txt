set(MODULE_NAME "CommonUtils")

O2_SETUP(NAME ${MODULE_NAME})

set(SRCS
  src/TreeStream.cxx
  src/TreeStreamRedirector.cxx
  src/RootChain.cxx
  src/CompStream.cxx
  src/ShmManager.cxx
)

Set(HEADERS
  include/${MODULE_NAME}/TreeStream.h
  include/${MODULE_NAME}/TreeStreamRedirector.h
  include/${MODULE_NAME}/RootChain.h
  include/${MODULE_NAME}/BoostSerializer.h
  include/${MODULE_NAME}/ShmManager.h
  include/${MODULE_NAME}/RngHelper.h
)

Set(LINKDEF src/CommonUtilsLinkDef.h)
Set(LIBRARY_NAME ${MODULE_NAME})
set(BUCKET_NAME common_utils_bucket)

O2_GENERATE_LIBRARY()

set(TEST_SRCS
  test/testTreeStream.cxx
  test/testBoostSerializer.cxx
  test/testCompStream.cxx
)

O2_GENERATE_TESTS(
  MODULE_LIBRARY_NAME ${LIBRARY_NAME}
  BUCKET_NAME ${BUCKET_NAME}
  TEST_SRCS ${TEST_SRCS}
)
