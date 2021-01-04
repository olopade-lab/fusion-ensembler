def flatten(data):
    result = data.copy()
    result[result > 1] = 1

    return result

def noop(data):
    return data
