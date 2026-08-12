# spark_parse() parses timeseries=FALSE calibration data

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["dim", "names", "numeric_sums", "numeric_na_counts"]
        }
      },
      "value": [
        {
          "type": "integer",
          "attributes": {},
          "value": [96, 19]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["calibrant", "fluorophore", "media", "concentration", "replicate", "well", "OD600", "OD700", "GFP 40", "GFP 50", "GFP 60", "GFP 70", "GFP 80", "GFP 90", "GFP 100", "GFP 110", "GFP 120", "row", "column"]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "replicate", "OD600", "OD700", "GFP 40", "GFP 50", "GFP 60", "GFP 70", "GFP 80", "GFP 90", "GFP 100", "GFP 110", "GFP 120", "column"]
            }
          },
          "value": [6021066735620238, 240, 17.54, 15.526, 15042, 80075, 303903, 447492, 293985, 342834, 366852, 384181, 395755, 624]
        },
        {
          "type": "integer",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "replicate", "OD600", "OD700", "GFP 40", "GFP 50", "GFP 60", "GFP 70", "GFP 80", "GFP 90", "GFP 100", "GFP 110", "GFP 120", "column"]
            }
          },
          "value": [0, 0, 0, 0, 0, 0, 0, 4, 12, 16, 20, 24, 28, 0]
        }
      ]
    }

# spark_parse() parses timeseries=TRUE sample data

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["dim", "names", "numeric_sums", "numeric_na_counts"]
        }
      },
      "value": [
        {
          "type": "integer",
          "attributes": {},
          "value": [4608, 20]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["strain", "host", "plasmid", "plasmid_2", "strain_2", "media", "sugar", "amino_acids", "inducer", "concentration", "init_ratio", "init_dilution", "well", "time", "OD600", "OD700", "GFP", "mCherry", "row", "column"]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "init_ratio", "init_dilution", "time", "OD600", "OD700", "GFP", "mCherry", "column"]
            }
          },
          "value": [0, 432, 61.383, 129937169.28, 343.887, 283.543, 24625314, 5820116, 29952]
        },
        {
          "type": "integer",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "init_ratio", "init_dilution", "time", "OD600", "OD700", "GFP", "mCherry", "column"]
            }
          },
          "value": [3600, 3744, 3744, 0, 3600, 3600, 3600, 3600, 0]
        }
      ]
    }

