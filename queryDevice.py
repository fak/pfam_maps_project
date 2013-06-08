
from django.db import connection

def custom_sql(query, params):
    cursor = connection.cursor()
    cursor.execute(query, params)
    data = cursor.fetchall()
    return data


