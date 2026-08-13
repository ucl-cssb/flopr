# generate_cfs() produces conversion factors from Tecan Spark calibration data

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
          "value": [12, 5]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["cf", "beta", "calibrant", "fluorophore", "measure"]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["cf", "beta"]
            }
          },
          "value": [0, 0]
        },
        {
          "type": "integer",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["cf", "beta"]
            }
          },
          "value": [0, 0]
        }
      ]
    }

# process_plate() handles multi-fluorophore + to_MEFL=TRUE (backlog 004, the README's own example)

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
          "value": [4608, 25]
        },
        {
          "type": "character",
          "attributes": {},
          "value": ["strain", "host", "plasmid", "plasmid_2", "strain_2", "media", "sugar", "amino_acids", "inducer", "concentration", "init_ratio", "init_dilution", "well", "time", "OD600", "OD700", "GFP", "mCherry", "row", "column", "normalised_OD700", "normalised_GFP", "normalised_mCherry", "calibrated_OD700", "calibrated_GFP"]
        },
        {
          "type": "double",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "init_ratio", "init_dilution", "time", "OD600", "OD700", "GFP", "mCherry", "column", "normalised_OD700", "normalised_GFP", "normalised_mCherry", "calibrated_OD700", "calibrated_GFP"]
            }
          },
          "value": [0, 432, 61.383, 129937169.28, 343.887, 283.543, 24625314, 5820116, 29952, 196.787, 11394270.493, 4885909.855, 169385799123.42001, 1047470878188852.2]
        },
        {
          "type": "integer",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["concentration", "init_ratio", "init_dilution", "time", "OD600", "OD700", "GFP", "mCherry", "column", "normalised_OD700", "normalised_GFP", "normalised_mCherry", "calibrated_OD700", "calibrated_GFP"]
            }
          },
          "value": [3600, 3744, 3744, 0, 3600, 3600, 3600, 3600, 0, 3600, 3600, 3600, 3600, 3600]
        }
      ]
    }

