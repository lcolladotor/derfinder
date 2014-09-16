## Run the tests if the system variable 'RUN_ALL_TESTS' is set to TRUE

flag <- as.logical(Sys.getenv('RUN_ALL_TESTS'))
if(!is.na(flag)) {
    if(flag == TRUE) {
        library('testthat')
        test_check('derfinder')
    }
} else {
    message('Set the system variable RUN_ALL_TESTS to TRUE to run all the tests.')
}
