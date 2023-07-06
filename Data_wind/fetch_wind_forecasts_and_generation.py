import requests
import io
import pandas as pd


# Function for data fetching
def get_data(url, headers, parameters, print_status=False ):
    api_response = requests.get(url, headers = headers, params = parameters)
    if print_status:
        print(api_response.status_code)

    return api_response.text


## -- Wind power generation hourly forecast - updated once a day -- 
dataset_ID_forecast = "246"
url_forecast = "https://api.fingrid.fi/v1/variable/"+dataset_ID_forecast+"/events/csv"

## -- Wind power hourly generation -- 
dataset_ID_generation = "75"
url_generation = "https://api.fingrid.fi/v1/variable/"+dataset_ID_generation+"/events/csv"

## -- Total production capacity used in the wind power forecast --
dataset_ID_capacity = "268"
url_capacity = "https://api.fingrid.fi/v1/variable/"+dataset_ID_capacity+"/events/csv"


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

    ## -- Wind power generation hourly forecast - updated once a day -- 
    wind_forecast = get_data(url_forecast, headers, parameters)

    # Save data to dataframe
    df_forecast = pd.read_csv(io.StringIO(wind_forecast))
    df_forecast.rename(columns={'value': 'forecast'}, inplace=True)
    #print(df_forecast)

    df_forecast['day'] = end_dates[d]
    df_forecast['hour'] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]



    ## -- Wind power hourly generation -- 
    wind_generation = get_data(url_generation, headers, parameters)

    # Save data to dataframe
    df_generation = pd.read_csv(io.StringIO(wind_generation))
    df_generation.rename(columns={'value': 'generation'}, inplace=True)
    #print(df_generation)



    ## -- Total production capacity used in the wind power forecast --
    wind_capacity = get_data(url_capacity, headers, parameters)
    
    # Save data to dataframe
    df_capacity = pd.read_csv(io.StringIO(wind_capacity))
    df_capacity.rename(columns={'value': 'capacity'}, inplace=True)
    #print(df_capacity)


    # # Join data frames and calculate wind factor
    df = df_forecast.set_index(['start_time', 'end_time']).join(df_capacity.set_index(['start_time', 'end_time']), lsuffix='_forecast', rsuffix = '_capacity')
    df['wind_factor_forecast'] = df['forecast']/df['capacity']
    #print(df)

    df = df.join(df_generation.set_index(['start_time', 'end_time']), lsuffix='_1', rsuffix = '_generation')
    df['wind_factor_generation'] = df['generation']/df['capacity']
    #print(df)

    df.set_index(['day', 'hour'], inplace=True)

    # Saving dataframe as csv
    df.to_csv("wind_data_"+dates[d]+'.csv')

