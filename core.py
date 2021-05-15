#core.py

import importlib
import re

class app:
    def __init__(self, config_path):
        self.metadata = None
        self.plugins = None
        self.handles = None
        self.workflow = None

        config = self.read_config(config_path)

        self.set_metadata(config["metadata"])
        self.set_plugins(config["plugins"])
        self.set_handles(config["handles"])
        self.set_workflow(config["workflow"])

    """
    This function reads the config file which is divided into four fields :
        - metadata: contains the general and shared informations for EMBLish
            as a couple key:value
        - plugins: contains the list of the plugins that will be used as a
            triplet plugin_key:plugin_name,plugin_package
        - handles: contains the list of files that will be used as inputs as
            a triplet handle_key:plugin,file_path
        - workflow: contains a hierarchical list of the different step to run
            as a couple plugin_key,handle_key
    """
    def read_config(self, config_path):
        config = {
            "metadata":[],
            "plugins":[],
            "handles":[],
            "workflow":[]
        }
        
        current_field = None
        current_line = None

        with open(config_path) as handle:
            for line in handle:
                current_line = line.rstrip("\n")

                #Comments and blank lines ignore
                if current_line == "": continue
                if current_line[0] == "#": continue

                #Filling the config object with the content of the configuration file
                if current_field:
                    assert current_field in ["metadata", "plugins", "handles", "workflow"]

                    #Checking if we are at the end of a field
                    if re.match(r"^(<\/)(\w+)(>)$", current_line):
                        current_field = None
                    else:
                        config[current_field].append(current_line)

                else:
                    assert re.match(r"^(<)(\w+)(>)$", current_line)
                    assert current_line[1:-1] in ["metadata", "plugins", "handles", "workflow"]

                    current_field = current_line[1:-1]

            handle.close()
        return config

    """
    This function converts the array containing the metadata into a dictionnary
    """
    def set_metadata(self, array):
        self.metadata = {element.split(":")[0]:element.split(":")[1] for element in array}

    """
    This function converts the array containing the plugins parameters into a dictionnary 
    with plugins to call with their key
    """
    def set_plugins(self, array):
        self.plugins = {element.split(":")[0]:importlib.import_module(element.split(":")[1].split(",")[0],element.split(":")[1].split(",")[1]).Plugin() for element in array}

    """
    This function converts the array containing the handles parameters into a dictionnary
    with handles to call with their key
    """
    def set_handles(self, array):
        self.handles = {element.split(":")[0]:self.plugins[element.split(":")[1].split(",")[0]].process(element.split(":")[1].split(",")[1]) for element in array}
    
    """
    """
    def set_workflow(self, array):
        temp = self.refactor_workflow(array)
        temp = self.merge_workflow(temp)
        self.workflow = self.convert_workflow_task(temp[0])

    """
    Convert the list element in triplet level,<plugin_key,handle_key>,[]
    """
    def refactor_workflow(self, array):

        for i in range(len(array)):
            element = array[i]
            level = len(re.search(r"^(-)+", element).group())

            array[i] = (level, element[level:], [])

        return array

    """
    Order the elements and create the hierarchical nodes
    """
    def merge_workflow(self, array):
        array.insert(0, (0,None,[]))

        max_level = max(array, key = lambda element: element[0])[0]

        for i in range(max_level, -1, -1):
            for j in range(len(array)):
                element = array[j]
    
                if element[0] == i:
                    for k in range(1, j+1):
                        if (
                            array[j-k] is not None and
                            array[j-k][0] == i-1
                        ):
                            array[j-k][2].append(element)
                            array[j] = None
                            break

            array = [element for element in array if element != None]
            
        return array

    """
    """
    def convert_workflow_task(self, task):
        if task[1]:
            return (
                self,
                task[1].split(",")[0],
                task[1].split(",")[1],
                [self.convert_workflow_task(sub_task) for sub_task in task[2]]
            )

        return [self.convert_workflow_task(sub_task) for sub_task in task[2]]
        
    """
    """
    def run(self):
        for app, key_plugin, *args in self.workflow:
            app.plugins[key_plugin].process(app, *args)