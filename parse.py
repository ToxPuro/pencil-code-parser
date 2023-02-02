
import math
import re
import os
import csv



class Parser:

    def __init__(self):
        self.static_variables = {}
        self.lines = {}
        self.parsed_files_for_static_variables = []
        self.parsed_modules = []
        self.parsed_subroutines = []
        self.loaded_files = []
        self.subroutine_order = 0
        self.used_files = []
        with open("make-output.txt", mode="r") as file:
            lines = file.readlines()
            for line in lines:
                search = re.search("([^\s]+\.f90)", line)
                if search is not None:
                    self.used_files.append(f"./pencil-code/src/{search.group(1)}")

    def find_module_files(self, module_name, rootdir):
        res = [file for file in self.used_files if self.contains_module(self.get_lines(file), module_name)]
        self.loaded_files.extend(res)
        self.loaded_files = list(set(self.loaded_files))
        return res

    def find_subroutine_files(self, subroutine_name, rootdir):
        return [x for x in self.used_files if self.contains_subroutine(self.get_lines(x), subroutine_name)]

    def parse_module(self, module_name, rootdir):
        if module_name not in self.parsed_modules:
            self.parsed_modules.append(module_name)
            file_paths = self.find_module_files(module_name, rootdir)
            for file_path in file_paths:
                self.parse_file_for_static_variables(file_path)


    def parse_line(self, line):
        ## remove comment at end of the line
        iter_index = len(line)-1;
        num_of_single_quotes = 0
        num_of_double_quotes = 0
        for iter_index in range(len(line)):
            if line[iter_index] == "'":
                num_of_single_quotes += 1
            if line[iter_index] == '"':
                num_of_double_quotes += 1
            if line[iter_index] == "!" and num_of_single_quotes%2 == 0 and num_of_double_quotes%2 == 0 and line[iter_index+1] != "=":
                line = line[0:iter_index]
                break                          
        return line.strip()

    def get_lines(self, filepath, start=0, end=math.inf, include_comments=False):
        if filepath not in self.lines.keys():
            lines = []
            read_lines = open(filepath, 'r').readlines()
            index = 0
            while index<len(read_lines):
                line = read_lines[index]
                if len(line)>0:
                    if line[0] != "!":
                        line = self.parse_line(line)
                    if len(line)>0:
                        if line[-1] == "&":
                            while line[-1] == "&":
                                index += 1
                                next_line = read_lines[index].strip()
                                if next_line[0] != "!" and len(next_line)>0:
                                    line = (line[:-1] + " " + self.parse_line(next_line)).strip()
                        lines.append((line, index))
                index += 1
            self.lines[filepath] = lines
        res = [x for x in self.lines[filepath] if x[1] >= start and x[1] <= end]
        if not include_comments:
            res = [x for x in res if x[0][0] != "!"] 
        return res

    def contains_module(self, lines, module_name):
        for (line, count) in lines:
            if line.strip() == f"module {module_name}":
                return True
        return False

    def contains_subroutine(self, lines, subroutine_name):
        for (line, count) in lines:
            if line.split("(")[0].strip() == f"subroutine {subroutine_name}":
                return True
        return False

    def get_used_modules(self,lines):
        modules = []
        for (line, count) in lines:
            if line.strip().split(" ")[0].strip() == "use":
                modules.append(line.strip().split(" ")[1].strip().replace(",",""))
        return modules

    def get_mem_access_or_function_calls(self,lines):
        res = []
        for(line, count ) in filter(lambda x: not self.is_variable_line(x),lines):
            line = line.strip()
            matches = re.findall("[^'=' '\/+-.*()<>]+\(.+\)", line)
            if len(matches) > 0:
                res.extend(matches)
        return res


    def parse_write_variable(self, line_segment):
        iter_index = len(line_segment)-1
        end_index = iter_index
        start_index = 0
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        while iter_index>0 and (line_segment[iter_index] in " ()" or num_of_left_brackets != num_of_right_brackets):
            elem = line_segment[iter_index]
            if elem == "(":
                num_of_left_brackets += 1
            elif elem == ")":
                num_of_right_brackets += 1
            iter_index -= 1
        end_index = iter_index+1

        while iter_index>0 and line_segment[iter_index] not in " *+-();":
            iter_index -= 1
        
        if iter_index == 0:
            res = line_segment[0:end_index]
        else:
            res = line_segment[iter_index+1:end_index]
        if "%" in res:
            end_index = 0
            while res[end_index] != "%":
                end_index += 1
            res = res[0:end_index]
        return res


    def get_writes_from_line(self,line,count):
        res = []
        index = 0
        start_index = 0
        num_of_single_quotes = 0
        num_of_double_quotes = 0
        num_of_left_brackets = 0
        num_of_right_brackets = 0
        while index<len(line)-1:
            index += 1
            elem = line[index]
            if elem == "'":
                num_of_single_quotes += 1
            elif elem == '"':
                num_of_double_quotes += 1
            elif elem == "=" and line[index-1] not in "<>=!" and line[index+1] not in "<>=!" and num_of_single_quotes%2==0 and num_of_double_quotes%2==0 and num_of_left_brackets == num_of_right_brackets:
                write = self.parse_write_variable(line[start_index:index])
                res.append((write, count, line))
                start_index = index
            elif elem == "(":
                num_of_left_brackets += 1
            elif elem == ")":
                num_of_right_brackets += 1
        return res

    def get_writes(self,lines):
        res = []
        for(line, count) in filter(lambda x: not self.is_variable_line(x),lines):
            res.extend(self.get_writes_from_line(line,count))
        return res



    def get_function_calls(self,lines, local_variables):
        function_calls = []
        for(line, count ) in filter(lambda x: not self.is_variable_line(x),lines):
            regex_matches = re.finditer("([^'=' '\/+-.*\(\)<>]+)\(", line)
            for match in regex_matches:
                function_call = match
                current_index = function_call.start(0)
                while(line[current_index]!="("):
                    current_index +=1
                function_name = line[function_call.start(0):current_index]
                if function_name not in self.static_variables and function_name not in local_variables:
                    #step through (
                    current_index +=1
                    number_of_right_brackets = 1
                    number_of_left_brackets = 0
                    parameter_list_start_index = current_index
                    while(number_of_right_brackets>number_of_left_brackets):
                        if current_index >= len(line):
                            print(line)
                        if line[current_index] == "(":
                            number_of_right_brackets += 1
                        elif line[current_index] == ")":
                            number_of_left_brackets += 1
                        current_index += 1
                        
                    parameter_list = line[parameter_list_start_index:current_index-1]
                    parameters = [x.strip() for x in parameter_list.split(",")]
                    function_calls.append({"function_name": function_name, "parameters": parameters})

                #    


        return function_calls



    def get_subroutine_line_start_and_end(self,filename, function_name):
        subroutine_line_start = 0
        subroutine_line_end = 0
        for (line, count) in self.get_lines(filename):
            if line.split("(")[0].strip() == f"subroutine {function_name}":
                subroutine_line_start = count

            # line.strip() == f"endsubroutine {function_name}" or line.strip() == f"end subroutine {function_name}"    
            if  (re.match(f"end.+{function_name}",line) or (line == "end" or line == "endsubroutine")) and subroutine_line_start>0:
                subroutine_line_end = count
                return (subroutine_line_start, subroutine_line_end)
        print(filename, function_name)
        return None

    def get_contains_line_num(self, filename):
        lines = self.get_lines(filename)
        for (line, count) in lines:
            if line.strip() == "contains":
                return count
            elif re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line):
                return count-1
        return len(lines)


    def is_variable_line(self, line_elem):
        parts = line_elem[0].split("::")
        if "!" in parts[0]:
            return False
        return len(parts)>1


    def get_variable_names_from_line(self,line_variables):
        res = []
        start_index=0
        end_index=0
        current_index=0
        parse_still = True
        while parse_still:
            current_elem = line_variables[current_index]
            if current_elem == "=":
            
                res.append(line_variables[start_index:current_index])
                parse_until_next_variable = True
                current_index +=1
                num_of_left_brackets = 0
                num_of_right_brackets = 0
                while parse_until_next_variable:
                    current_elem = line_variables[current_index]
                    not_inside_brackets = num_of_left_brackets == num_of_right_brackets
                    if current_elem == "!":
                        parse_until_next_variable = False
                        parse_still = False
                    if current_elem == "," and not_inside_brackets:
                        start_index = current_index+1
                        parse_until_next_variable=False
                    if current_elem == "(":
                        num_of_left_brackets += 1
                    if current_elem == ")":
                        num_of_right_brackets += 1
                    if parse_until_next_variable:
                        current_index +=1
                        if current_index >= len(line_variables):
                            start_index = current_index+1
                            parse_until_next_variable = False
                            parse_still = False
            elif current_elem == ",":
                res.append(line_variables[start_index:current_index])
                start_index = current_index + 1
            elif current_elem == "!":
                parse_still = False
                res.append(line_variables[start_index:current_index+1])
            current_index += 1
            if current_index >= len(line_variables):
                parse_still= False
                if current_index > start_index:
                    res.append(line_variables[start_index:current_index+1])
        return [x.replace("&","").replace("!","").replace("(:)","").strip() for x in res if x.replace("&","").replace("!","").strip() != ""]

    def get_variables(self, lines, variables, filename):
        in_type = False
        for (line, count) in filter(self.is_variable_line, lines):
            parts = line.split("::")
            is_variable = True
            start =  parts[0].strip()
            type = start.split(",")[0].strip()
            if type == "public":
                is_variable = False
            if line.split(" ")[0] == "type" and len(line.split("::")) == 1:
                in_type = True
            if line.split(" ")[0] == "endtype":
                in_type = False
            if is_variable and not in_type:
                allocatable = "allocatable" in [x.strip() for x in start.split(",")]
                public = "public" in [x.strip() for x in start.split(",")]
                dimension = re.search("dimension\s*\((.+?)\)", start)
                if dimension is None:
                    dimension = []
                else:
                    dimension = [x.strip() for x in dimension.group(1).split(",")]
                line_variables = parts[1].strip()
                variable_names = self.get_variable_names_from_line(line_variables)
                for variable_name in variable_names:
                    variables[variable_name] = {"type": type, "dims": dimension, "allocatable": allocatable, "origin": filename, "public": public, "threadprivate": False}
        return variables


    def parse_file_for_static_variables(self, filepath):
        if filepath not in self.parsed_files_for_static_variables:
            self.parsed_files_for_static_variables.append(filepath)
            modules = self.get_always_used_modules(filepath)
            for module in modules:
                self.parse_module(module, "./pencil-code/src")
            self.load_static_variables(filepath)

    def parse_subroutine_for_static_variables(self,subroutine_name, filepath):
        self.parse_file_for_static_variables(filepath)
        for module in self.get_subroutine_modules(filepath, subroutine_name):
            self.parse_module(module, "./pencil-code/src")

    def add_threadprivate_declarations_in_file(self,filename):
        lines = self.get_lines(filename, include_comments=True)
        res = []
        for (line,count) in lines:
            is_threadprivate_declaration = "!$omp" in line and "threadprivate" in line
            if is_threadprivate_declaration:
                variable_names = [variable.strip() for variable in re.search("threadprivate\((.+)\)",line).group(1).split(",")]
                for variable in variable_names:
                    if variable in self.static_variables:
                        self.static_variables[variable]["threadprivate"] = True
                        res.append(variable)
        return res

    def add_public_declarations_in_file(self,filename,lines):
        for (line,count) in lines:
            parts = line.split("::")
            is_public_declaration = "public" == parts[0].split(",")[0].strip() and len(parts) == 2
            if is_public_declaration:
                variable_names = self.get_variable_names_from_line(parts[1])
                for variable in variable_names:
                    if variable in self.static_variables:
                        self.static_variables[variable]["public"] = True 
            match = re.search("include\s+(.+\.h)",line)
            if match:
                header_filename = match.group(1).replace("'","").replace('"',"")
                directory = re.search('(.+)\/',filename).group(1)
                header_filepath = f"{directory}/{header_filename}"
                with open(header_filepath,"r") as file:
                    lines = file.readlines()
                    for line in lines:
                        parts = line.split("::")
                        is_public_declaration = "public" == parts[0].split(",")[0].strip()
                        if is_public_declaration:
                            variable_names = self.get_variable_names_from_line(parts[1])
                            for variable in variable_names:
                                if variable in self.static_variables:
                                    self.static_variables[variable]["public"] = True 


    def load_static_variables(self, filename):
        static_variables_end = self.get_contains_line_num(filename)
        self.get_variables(self.get_lines(filename, 1, static_variables_end), self.static_variables, filename)
        self.add_public_declarations_in_file(filename,self.get_lines(filename, 1, static_variables_end))
        self.add_threadprivate_declarations_in_file(filename)

    def get_subroutine_variables(self, filename, subroutine_name):
        subroutine_line_start,subroutine_line_end = self.get_subroutine_line_start_and_end(filename, subroutine_name)
        return self.get_variables(get_lines(filename, subroutine_line_start,subroutine_line_end),filename)

    def get_always_used_modules(self, filename):
        module_declaration_end = self.get_contains_line_num(filename)
        return self.get_used_modules(self.get_lines(filename, 1, module_declaration_end))

    def get_subroutine_modules(self, filename, subroutine_name):
        subroutine_line_start, subroutine_line_end = self.get_subroutine_line_start_and_end(filename, subroutine_name)
        return self.get_used_modules(self.get_lines(filename, subroutine_line_start, subroutine_line_end))

    def get_subroutine_lines(self, filename, subroutine_name):
        subroutine_line_start, subroutine_line_end = self.get_subroutine_line_start_and_end(filename, subroutine_name)
        return self.get_lines(filename, subroutine_line_start, subroutine_line_end)

    def parse_subroutine_all_files(self, subroutine_name, rootdir, call_trace, layer_depth=math.inf):
        self.subroutine_order += 1
        static_writes = []
        if subroutine_name not in self.parsed_subroutines:
            self.parsed_subroutines.append(subroutine_name)
            file_paths = self.find_subroutine_files(subroutine_name, rootdir)
            for file_path in file_paths:
                static_writes.extend(self.parse_subroutine_in_file(file_path, subroutine_name, layer_depth, call_trace))
        return static_writes


    def parse_subroutine_in_file(self, filename, subroutine_name, layer_depth=math.inf, call_trace=""):

        if layer_depth < 0:
            return []
        subroutine_line_start, subroutine_line_end = self.get_subroutine_line_start_and_end(filename, subroutine_name)
        lines = self.get_lines(filename, subroutine_line_start, subroutine_line_end)
        modules = self.get_used_modules(lines)
        self.parse_file_for_static_variables(filename)
        for module in self.get_subroutine_modules(filename, subroutine_name):
            self.parse_module(module, "./pencil-code/src")

        local_variables = {}
        local_variables = self.get_variables(lines, local_variables,filename)
        function_calls = self.get_function_calls(lines, local_variables)
        
        writes = self.get_writes(lines)
        if self.subroutine_order == 0:
            call_trace = subroutine_name
        else:
            call_trace = f"{call_trace} -> {subroutine_name}"
        static_writes = [{"variable": x[0],"line_num": x[1], "filename": filename, "call_trace":call_trace, "line": x[2]} for x in writes if (x[0] not in local_variables and x[0] in self.static_variables)]
        for function_call in function_calls:
            static_writes.extend(self.parse_subroutine_all_files(function_call["function_name"], "./pencil-code/src", call_trace, layer_depth-1))
        return static_writes

    def generate_save_array_store(self,store_variable):
        res = f"{store_variable}_generated_array(imn"
        if len(self.static_variables[store_variable]['dims']) == 0:
            res += ",1"
        else:
            for i in range(len(self.static_variables[store_variable]['dims'])):
                res += ",:"
        res += f") = {store_variable}\n"
        return res

    def generate_read_from_save_array(self,store_variable):
        res = f"{store_variable} = {store_variable}_generated_array(imn"
        if len(self.static_variables[store_variable]['dims']) == 0:
            res += ",1"
        else:
            for i in range(len(self.static_variables[store_variable]['dims'])):
                res += ",:"
        res +=")\n"
        return res

    def generate_allocation_for_save_array(self,store_variable):
        res = f"{self.static_variables[store_variable]['type']}, dimension ("
        if self.static_variables[store_variable]["allocatable"]:
            res += ":"
        else:
            res += "nx*ny"
        if len(self.static_variables[store_variable]['dims']) == 0:
            if self.static_variables[store_variable]["allocatable"]:
                res += ",:"
            else:
                res += ",1"
        else:
            for dimension in self.static_variables[store_variable]["dims"]:
                res += f",{dimension}"
        res += ")"
        if self.static_variables[store_variable]["allocatable"]:
            res += ", allocatable"
        res += f" :: {store_variable}_generated_array\n"
        return res

    def save_static_variables(self):
        with open("static_variables.csv","w",newline='') as csvfile:
            writer = csv.writer(csvfile)
            for variable in self.static_variables.keys():
                writer.writerow([variable,self.static_variables[variable]["type"],self.static_variables[variable]["dims"],self.static_variables[variable]["allocatable"],self.static_variables[variable]["origin"],self.static_variables[variable]["public"],self.static_variables[variable]["threadprivate"]])
        
    def read_static_variables(self):
        with open("static_variables.csv","r",newline='') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                self.static_variables[row[0]] = {"type": row[1], "dims": [dim for dim in row[2].replace("'","").strip('][').split(', ') if dim != ""], "allocatable": (row[3].lower() in ("yes", "true", "t", "1")), "origin": row[4], "public": (row[5].lower() in ("yes", "true", "t", "1")), "threadprivate": (row[6].lower() in ("yes", "true", "t", "1"))}

    def save_static_writes(self,static_writes):
        keys = static_writes[0].keys()
        with open("writes.csv","w",newline='') as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(static_writes)

    def read_static_writes(self):
        with open("writes.csv", mode="r") as infile:
            return [dict for dict in csv.DictReader(infile)]


    def make_variables_threadprivate(self,variables):
        for file in self.used_files:
            already_threadprivate = self.add_threadprivate_declarations_in_file(file)
            res_contents = []
            contents = self.get_lines(file,include_comments=True)
            in_module_declaration = True
            in_type = False
            for(line,count) in contents:
                res_contents.append(line)
                if line.split(" ")[0] == "type":
                    in_type = True
                if line.split(" ")[0] == "endtype":
                    in_type = False
                if in_module_declaration and not in_type:
                    if re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line) or line.strip() == "contains":
                        in_module_declaration = False
                    parts = line.split("::")
                    if len(parts) > 1:
                        line_variables = parts[1].strip()
                        variable_names = self.get_variable_names_from_line(line_variables)
                        res_contents.extend([f"!$omp threadprivate({variable})" for variable in variable_names if variable in variables and variable not in already_threadprivate])
            
            with open(f"./out/{file}","w") as f:
                f.write("\n".join(res_contents))
                    
    def add_atomic_declarations(self,variables):
        files = self.used_files
        for file in files:
            res_contents = []
            contents = self.get_lines(file,include_comments=True)
            in_module_declaration = True
            in_type = False
            for (line,count) in contents:
                res_contents.append(line)
                writes = self.get_writes_from_line(line,count)

                #handle conditional write and where writes
                if "if" in line or line.split("(")[0].strip() == "where":
                    for variable in variables:
                        if variable in [write[0] for write in writes]:
                            last_line = res_contents[-1]
                            res_contents[-1] = "!$omp critical\n"
                            res_contents.append(last_line)
                            res_contents.append("!$omp end critical\n")
                else:
                    for variable in variables:
                        if variable in [write[0] for write in writes]:
                            last_line = res_contents[-1]
                            res_contents[-1] = "!$omp critical\n"
                            res_contents.append(last_line)
                            res_contents.append("!$omp end critical\n")
                            #If one one's the more performant omp atomic
                            #last_line = res_contents[-1]
                            #res_contents[-1] = "!$omp atomic\n"
                            #res_contents.append(last_line)

            with open(f"./out/{file}","w") as f:
                f.write("\n".join(res_contents))

    def make_variables_public(self,variables):
        files = list(set([self.static_variables[variable]["origin"] for variable in variables ]))
        for file in files:
            res_contents = []
            contents = self.get_lines(file,include_comments=True)
            in_module_declaration = True
            in_type = False
            for (line,count) in contents:
                res_contents.append(line)
                if line.split(" ")[0] == "type":
                    in_type = True
                if line.split(" ")[0] == "endtype":
                    in_type = False
                if in_module_declaration and not in_type:
                    if re.match("function\s+.+\(.+\)",line) or re.match("subroutine\s+.+\(.+\)",line) or line.strip() == "contains":
                        in_module_declaration = False
                    parts = line.split("::")
                    if len(parts) > 1:
                        line_variables = parts[1].strip()
                        variable_names = self.get_variable_names_from_line(line_variables)
                        res_contents.extend([f"public :: {variable}" for variable in variable_names if variable in variables and self.static_variables[variable]["origin"] == file and not self.static_variables[variable]["public"]])

            with open(f"./out/{file}","w") as f:
                f.write("\n".join(res_contents))

        
        
    def generate_commands(self, filename, save_variables):
        res_contents = []
        with open(f"{filename}.f90", "r") as f:
            contents = f.readlines()
        for index in range(len(contents)):
            res_contents.append(contents[index])
            if re.match("!\$parser-command:save-global-state.?", contents[index]):
                save_string = "".join([self.generate_save_array_store(static_variable) for static_variable in save_variables])
                res_contents.append(save_string)
            elif re.match("!\$parser-command:load-global-state.?", contents[index]):
                save_string = "".join([self.generate_read_from_save_array(static_variable) for static_variable in save_variables])
                res_contents.append(save_string)
            elif re.match("!\$parser-command:allocate-global-state-arrays.?", contents[index]):
                save_string = "".join([self.generate_allocation_for_save_array(static_variable) for static_variable in save_variables])
                res_contents.append(save_string)
            elif re.match(".+!\$parser-command:generate-firstprivate-pragma-for-static-variables.?",contents[index]):

                #remove the parser command for the omp pragma to be valid
                res_contents[-1] = res_contents[-1].replace("!$parser-command:generate-firstprivate-pragma-for-static-variables","firstprivate(&")
                save_string = ",&\n".join([f"!$omp {variable}" for variable in save_variables])
                save_string = f"{save_string})\n"
                res_contents.append(save_string)
                


        with open(f"./out/{filename}.f90", "w") as f:
            contents = "".join(res_contents)
            f.write(contents)


def main():
    filename = "./pencil-code/src/equ.f90"
    subroutime_name = "rhs_cpu"
    parser = Parser()
    #parser.load_static_variables("./pencil-code/src/cdata.f90")
   # static_writes = parser.parse_subroutine_in_file("./pencil-code/src/entropy.f90","calc_heat_cool",0)
   # parser.save_static_writes(static_writes)
    #for (line,count) in lines:
    #    print(line)
    #parser.parse_file_for_static_variables("./pencil-code/src/hydro.f90")

   # lines = parser.get_lines("./pencil-code/src/geometrical_types.f90",include_comments=True)
   # for (line,count) in lines:
   #     print(line)

    #static_writes = parser.parse_subroutine_in_file(filename, subroutime_name)
    #parser.save_static_variables()
    #parser.save_static_writes(static_writes)
    static_writes = parser.read_static_writes()
    variables = list(set([x["variable"] for x in static_writes]))
    variables = list(set(variables) - set(variables[:len(variables)//2]))
    variables = list(set(variables) - set(["fname","fnamer","fnamex","fnamey","fnamez","fnamexy","fnamexz","fnamerz","itype_name","dt1_max_array","dt1_max"]))
    variables = variables[:len(variables)//2]
    for variable in variables:
        print(variable)
    #parser.make_variables_threadprivate(variables)
    #print(f"Num of static variables written with full depth {len(variables)}")
    #parser.add_atomic_declarations(["fname","fnamer","fnamex","fnamey","fnamez","fnamexy","fnamexz","fnamerz","itypename"])

    #filename = "./pencil-code/src/entropy.f90"
    #subroutine_name = "calc_heat_cool"
    #subroutine_line_start, subroutine_line_end = parser.get_subroutine_line_start_and_end(filename, subroutine_name)
    #lines = parser.get_lines(filename, subroutine_line_start, subroutine_line_end)
   # parser.parse_file_for_static_variables(filename)
   # for module in parser.get_subroutine_modules(filename, subroutine_name):
   #     parser.parse_module(module, "./pencil-code/src")

    #local_variables = {}
    #local_variables = parser.get_variables(lines, local_variables,filename)
    #for variable in local_variables:
    #    print(variable)
if __name__ == "__main__":
    main()