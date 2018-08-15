import sys
import zlib
import gzip
import zipfile
import os
import shutil


#will create file if it doesn't exist
z = zipfile.ZipFile("new_file.zip", "a", zipfile.ZIP_DEFLATED)
z.write(sys.argv[1])
z.close()