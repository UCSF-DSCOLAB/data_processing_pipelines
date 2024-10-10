class TestingFunctions {


    def static get_result_files(String testDirAbsolute) {
            def finding_doublets = "/finding_doublets"
            def automated_processing = "/automated_processing"
            def freemuxlet = "/freemuxlet"
            def test_directories = [finding_doublets, automated_processing, freemuxlet]

        // Setup assertion directories
            def dm1_scg1 = testDirAbsolute + "/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG1"
            def dm1_scg2 = testDirAbsolute + "/data/single_cell_GEX/processed/TEST-POOL-DM1-SCG2"
            def dm2_scg1 = testDirAbsolute + "/data/single_cell_GEX/processed/TEST-POOL-DM2-SCG1"
            def test_sample_paths = [dm1_scg1, dm1_scg2, dm2_scg1]
            def resultFilesList = []
            for (sample in test_sample_paths) {
                for (test_dir in test_directories) {
                    def dir_to_test = sample + test_dir
                    def dirFile = new File(dir_to_test)
                    def resultFiles = dirFile.listFiles().collect { it.name }
                    resultFilesList.addAll(resultFiles)
                }
            }

            // Filter out the .pdf files
            def filteredPaths = resultFilesList.findAll { filePath -> !filePath.endsWith('.pdf')}
            def sorted = filteredPaths.sort()

            return sorted
    }
}

