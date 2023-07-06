import requests
import pandas as pd
from bs4 import BeautifulSoup

# Function for data fetching
def get_data(url, headers, parameters, print_status=False ):
    response = requests.get(url, headers = headers, params = parameters)
    if print_status:
        print(response.status_code)

    soup = BeautifulSoup(response.content, 'xml')
    return soup


# Unpacking day-ahead prices
def unpack_DA_prices(soup):
    prices = {'day':[], 'hour':[], 'price':[]}

    for n, tag in enumerate(soup.find_all('Period')):
        day = tag.find('end').get_text()[0:10]
        for i, timeseries in enumerate(tag.find_all('Point')):
            hour = timeseries.find('position').get_text()
            price = timeseries.find('price.amount').get_text()
            prices['day'].append(day)
            prices['hour'].append(hour)
            prices['price'].append(price)
    return pd.DataFrame(prices)


####################################################################################################
# DAY-AHEAD PRICES
## -- Setting up REST GET--
end_point = "https://web-api.tp.entsoe.eu/api" 

# Specify dataset
documentType = "A44"

# Domain Finland
in_domain = "10YFI-1--------U"
out_domain = "10YFI-1--------U"

# Security token
token = "ADD SECURITY TOKEN HERE"

# Data needed
dates = ["20220301", "20220302", "20220308", "20220309", "20220315", "20220316",
        "20230221","20230222", "20230228", "20230301", "20230307"] 


# Fetch data for each date and save it in CSV
for d in dates:

    # Note: these are in CET (but return values have UTC printed out which correspond to these CET times). Format: yyyyMMddHHmm
    start_time = d+"0000"
    end_time = d+"2300"


    parameters = {"securityToken": token, "documentType": documentType, "in_Domain": in_domain, "out_Domain": out_domain, "periodStart": start_time, "periodEnd": end_time}
    headers={}



    # Fetch data
    soup = get_data(end_point, headers, parameters)

    # Print XML
    #print(soup.prettify())

    # Unpack XML to dataframe
    df = unpack_DA_prices(soup)


    # Save data to csv
    df.to_csv('DA_prices_'+d+'.csv', index=False)

####################################################################################################