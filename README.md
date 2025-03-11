# TIFF to COG Generation Tool

## Requirements

1. Python Requirements:
   - Python 3.7 or higher
   - GDAL library (python-gdal)
   - rasterio
   - watchdog
   - configparser

2. System Requirements:
   - Linux/Unix-based system
   - Sufficient disk space for TIFF and COG files
   - Read/Write permissions for:
     - Watch directories
     - Output directory
     - Logs directory

## Setup

1. Install system dependencies:
```bash
sudo apt-get update
sudo apt-get install gdal-bin python3-gdal
```

2. Install Python dependencies:
```bash
pip install rasterio watchdog configparser
```

3. Configure the config.ini file:
   - Set appropriate paths for WATCH_DIR
   - Set OUTPUT_DIR for COG files
   - Ensure LOGS directory exists
   - Configure appropriate satellite and processing codes

4. Ensure all directories specified in config.ini exist and have proper permissions:
```bash
mkdir -p /path/to/watch/dir1
mkdir -p /path/to/output/dir
mkdir -p /usr/local/logs/CogGen
chmod 755 /usr/local/logs/CogGen
```

## Note
- Ensure all specified directories in config.ini are accessible and have proper read/write permissions
- The watch directories should contain TIFF files matching the specified SAT_CODES and LPP_CODES
- Log directory must be writable by the process