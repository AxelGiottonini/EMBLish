#core.py

import importlib
import re

"""
"""
class App:
    """
        Description:

        Arguments:

        Output:

        Note:
    """
    def __init__(self, config_path):

        self.metadata = None
        self.plugins = None
        self.handles = None
        self.workflow = None

        self.current_sequence = None

        config = self.read_config(config_path)

        self.set_metadata(config["metadata"])
        self.set_plugins(config["plugins"])
        self.set_handles(config["handles"])
        self.set_workflow(config["workflow"])

        self.all_plugins_required_metadata_check()

    """
        Description:

        Arguments:

        Output:

        Note:
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
        Description:
            - converts the array containing the metadata into a dictionnary by splitting
            each string containing the metadata key and metadata value into an item with 
            the key and value.
        Arguments:
            - array: list of strings
        Output:
            - dictionnary of strings
        Note:
    """
    def set_metadata(self, array:list):
        
        def convert(value):
            temp = value.split(",")
            if len(temp) > 1:
                if temp[1] == "int":
                    return int(temp[0])
            return value

        self.metadata = {element.split(":")[0]:convert(element.split(":")[1]) for element in array}

    """
        Description:
            - converts the array containing the plugins parameters (name, package) into a
            dictionnary by splitting each string containing the plugin key, the plugin name
            and the plugin package into an item with the key and the callable plugin.
            Store the result in the self.metadata variable.
        Arguments:
            - array: list of strings
        Output:
            - dictionnary of plugin objects
        Note:
    """
    def set_plugins(self, array:list):
        self.plugins = {element.split(":")[0]:importlib.import_module(element.split(":")[1].split(",")[0],element.split(":")[1].split(",")[1]).Plugin() for element in array}

    """
        Description:
            - converts the array containing the handles parameters into a dictionnary by 
            splitting each string containing the handle key, the handle converter and the
            file path into an item with the key and the converted as a data frame handle.
            Store the result in the self.plugin variable.
        Arguments:
            - array: list of strings
        Output:
            - dictionnary of data frames
        Note:
    """
    def set_handles(self, array):
        self.handles = {element.split(":")[0]:self.plugins[element.split(":")[1].split(",")[0]].process(element.split(":")[1].split(",")[1]) for element in array}
    
    """
        Description:
            - converts the array containing the workflow into a recursive automaton where
            task are described by a tuple containing the required plugin, the handle where
            the data is found and a list of elements to call.
            Store the result in the self.handles variable.
        Arguments:
            - array: list of strings
        Output (assigned):
            - array: list of tuples (recursive)
        Note:
    """
    def set_workflow(self, array):
        temp = self.refactor_workflow(array)
        temp = self.merge_workflow(temp)
        self.workflow = self.convert_workflow_task(temp[0])

    """
        Description:
            - converts strings into a tuple containing the level of the task, the rest of the 
            string and an empty array.
        Arguments:
            - array list of strings
        Output:
            - array: list of tuples
        Note:
    """
    def refactor_workflow(self, array):

        for i in range(len(array)):
            element = array[i]
            level = len(re.search(r"^(-)+", element).group())

            array[i] = (level, element[level:], [])

        return array

    """
        Description:
            - place the tasks in their parent (level-1) tasks array 
        Arguments:
            - array: list of tuples
        Output:
            -  array: list of tuples (recursive)
        Note:
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
        Description:
            - recursively converts the tuples containing the task string into a tuple containing 
            a reference to the application, the plugin key, the handle key and an array of 
            subtasks.
        Arguments:
            - array: list of tuples (recursive)
        Output:
            - array: list of tuples (recusrive)
        Note:
    """
    def convert_workflow_task(self, task):
        if task[1]:
            return (
                self,
                task[1].split(",")[0],
                task[1].split(",")[2] if len(task[1].split(",")) > 2 else "default",
                task[1].split(",")[1],
                [self.convert_workflow_task(sub_task) for sub_task in task[2]]
            )

        return [self.convert_workflow_task(sub_task) for sub_task in task[2]]
        
    """
        Description:

        Arguments:

        Output:

        Note:
    """
    def all_plugins_required_metadata_check(self):
        for key, plugin in self.plugins.items():
            if not plugin.required_metadata_check(self):
                raise InvalidConfigurationError(f"{key} plugin could not find required metadata")

    """
        Description:

        Arguments:

        Output:

        Note:
    """
    def run(self):
        for app, key_plugin, *args in self.workflow:
            app.plugins[key_plugin].process(app, *args)

"""
"""
class InvalidConfigurationError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f"InvalidConfigurationError, {self.message}"
        else:
            return "InvalidConfigurationError has been raised"