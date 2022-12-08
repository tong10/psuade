add_library(svmlight ${LIBRARY_TYPE} svm_learn_main.c svm_learn.c svm_common.c svm_classify.c svm_hideo.c)  

install (TARGETS svmlight
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

