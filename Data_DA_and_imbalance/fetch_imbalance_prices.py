import requests
import io
import pandas as pd
import numpy as np


# Function for data fetching
def get_data(url, headers, parameters, print_status=False ):
    api_response = requests.get(url, headers = headers, params = parameters)
    if print_status:
        print(api_response.status_code)

    return api_response.text


## -- Balancing Capacity Market (mFRR) hourly market price, down-regulating and up-regulating prices -- 
dataset_ID_imbalance = "319"
url_imbalance = "https://api.fingrid.fi/v1/variable/"+dataset_ID_imbalance+"/events/csv"

## -- Down-regulating balancing market prices
dataset_ID_down = "106"
url_down = "https://api.fingrid.fi/v1/variable/"+dataset_ID_down+"/events/csv"

## -- Up-regulating balancing market prices
dataset_ID_up = "244"
url_up = "https://api.fingrid.fi/v1/variable/"+dataset_ID_up+"/events/csv"


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


    # Imbalance prices 
    data = get_data(url_imbalance, headers, parameters)

    imbalance_prices = pd.read_csv(io.StringIO(data))
    imbalance_prices['day'] = end_dates[d]
    imbalance_prices['hour'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]


    # Down regulation prices 
    down_data = get_data(url_down, headers, parameters)
    down = pd.read_csv(io.StringIO(down_data))


    # Up regulation prices 
    up_data = get_data(url_up, headers, parameters)
    up = pd.read_csv(io.StringIO(up_data))
    

    df = imbalance_prices.set_index(['start_time', 'end_time']).join(down.set_index(['start_time', 'end_time']), lsuffix='', rsuffix = '_down')
    df = df.join(up.set_index(['start_time', 'end_time']), lsuffix='', rsuffix = '_up')
    df.columns = ["imbalance_price", "day", "hour", "price_down", "price_up"]
    #print(df)
    df = df.set_index(['day', 'hour'])

    # Get DA prices
    DA_prices = pd.read_csv('DA_prices_'+dates[d]+'.csv')

    df = df.join(DA_prices.set_index(['day', 'hour']))
    df.columns =  ["imbalance_price",  "price_down", "price_up", 'DA_price']

    # Calculate imbalance premium and deduce system imbalance
    df['difference'] = round(df['DA_price'] -  df['imbalance_price'], 2)
    df['system_imbalance'] = np.where(df['difference'] == 0, "-", "Regulation")
    #print(df)
    
    df.to_csv('imbalance_prices_'+dates[d]+'.csv', index=True)



