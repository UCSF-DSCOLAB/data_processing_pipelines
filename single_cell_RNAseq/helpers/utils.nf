def extractFileName(String path) {
    def filename = path.split('/').last() // Split by '/' and get the last part which is the filename
    return filename.split("\\.")[0] // Split the filename on dot and return the first part
}
