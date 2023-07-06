import requests
import io
import pandas as pd


# Function for data fetching
def get_data(url, headers, parameters, print_status=False ):
    api_response = requests.get(url, headers = headers, params = parameters)
    if print_status:
        print(api_response.status_code)

    return api_response.text


## -- Up-regulating FCR-D hourly market prices -- 
# Dataset ID 
dataset_ID = "318"
url = "https://api.fingrid.fi/v1/variable/"+dataset_ID+"/events/csv"



# Access key
headers = {"x-api-key": "ADD ACCESS KEY HERE"}


# Data needed
dates = ["20220301", "20220302", "20220308", "20220309", "20220315", "20220316",
        "20230221","20230222", "20230228", "20230301", "20230307"] 


# These dates and times are in EET, which correspond to CET 00:00-23:00
start_dates = ["2022-02-28", "2022-03-01", "2022-03-07", "2022-03-08", "2022-03-14", "2022-03-15",
        "2023-02-20","2023-02-21", "2023-02-27", "2023-02-28", "2023-03-06"]

end_dates = ["2022-03-01", "2022-03-02", "2022-03-08", "2022-03-09", "2022-03-15", "2022-03-16",
        "2023-02-21","2023-02-22", "2023-02-28", "2023-03-01", "2023-03-07"]

CET_midnight_in_UTC = [23, 23, 23, 23, 23, 23,
                       23, 23, 23, 23, 23]


# Fetch data for each date and save it in CSV
for d in range(0,len(dates)):

    # Note: these are in EET. Format YYYY-MM-ddTHH:mm:ssZ
    start_time = start_dates[d]+"T"+str(CET_midnight_in_UTC[d])+":00:00Z"
    end_time = end_dates[d]+"T"+str(CET_midnight_in_UTC[d]-1)+":59:00Z"
    parameters = {"start_time": start_time, "end_time": end_time}

    # Fetch data using function get_data
    prices = get_data(url, headers, parameters)

    # Save data to dataframe
    df = pd.read_csv(io.StringIO(prices))
    df.columns = ["start_time_UTC","end_time_UTC","price"]
    df['day'] = df['end_time_UTC'].str.slice(0, 10)
    df['hour'] = df['end_time_UTC'].str.slice(11, 13).astype('int')+1


    df.to_csv('FCR_up_'+dates[d]+'.csv', index=False)





