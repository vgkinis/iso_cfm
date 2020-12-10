import json


def modify_json(json_filepath, **kwargs):
    """

    """
    with open(json_filepath, "r") as f:
        json_string = f.read()
        c = json.loads(json_string)

    c = json.loads(json_string)
    if len(kwargs.keys()) != 0:
        for i in kwargs.keys():
            if i in c.keys():
                save_json = True
                c[i] = kwargs[i]
                print("changing json parameter %s, new value %s" %(i, c[i]))
            else: save_json = False
    if save_json == True:
        with open(json_filepath, "w") as f:
            json.dump(c, f, indent = 4, sort_keys = True)
    return
