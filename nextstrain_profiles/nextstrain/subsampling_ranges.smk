import datetime

# Set subsampling max date to today.
today = datetime.date.today()

# Set the earliest date to roughly 4 months ago (18 weeks).
early_late_cutoff = today - datetime.timedelta(weeks=18)
recent_late_cutoff = today - datetime.timedelta(weeks=4)

for build in config["subsampling"]:
    for scheme in config["subsampling"][build]:
        if "_recent" in scheme:
            config["subsampling"][build][scheme]["min_date"] = f"--min-date {recent_late_cutoff.strftime('%Y-%m-%d')}"
        if "_early" in scheme:
            config["subsampling"][build][scheme]["max_date"] = f"--max-date {early_late_cutoff.strftime('%Y-%m-%d')}"
        if "_late" in scheme:
            config["subsampling"][build][scheme]["min_date"] = f"--min-date {early_late_cutoff.strftime('%Y-%m-%d')}"
            config["subsampling"][build][scheme]["max_date"] = f"--max-date {recent_late_cutoff.strftime('%Y-%m-%d')}"
