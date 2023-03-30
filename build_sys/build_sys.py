
## This is the set of targets that need to be built 
targets = {'a': {'dependencies': ['b', 'c'], 'build_time': 10},
           'b': {'dependencies': ['d', 'c'], 'build_time': 1},
           'c': {'dependencies': None, 'build_time': 5},
           'd': {'dependencies': ['e', 'f', 'g'], 'build_time': 2},
           'e': {'dependencies': None, 'build_time': 9},
           'f': {'dependencies': None, 'build_time': 8},
           'g': {'dependencies': None, 'build_time': 7}}

# This example describes the following build:
# Target 'a' depends on Target 'b' and Target 'c'
# Target 'b' depends on Target 'd' and Target 'c'
# Target 'c' has no dependencies
# Target 'd' depends on Target 'e', Target 'f', Target 'g'
# Target 'e' has no dependencies
# Target 'f' has no dependencies
# Target 'g' has no dependencies


## Class object for single build target
class BuildTarget:
    ### You can modify stuff in here

    # Constructor for BuildTarget
    #
    #   return -> BuildTarget (instance of this class)
    #
    #   @param name The name of the target
    #   @param build_time The amount of time that target needs to build
    def __init__(self, name, build_time):
        # You can modify stuff in here
        self._name = name
        self._build_time = build_time

    # Simple function to build/display the output of building the target
    #
    #   return -> None
    #
    def build(self):
        ### You can modify stuff in here

        #### Do not modify these lines ###
        print("Built Target: '" + str(self._name) + "' Took: " + str(self._build_time) + "s")

    # Feel free to add methods to this class if needed


## Class object for the logic of building all targets
class BuildScript:
    ### You can modify stuff in here

    # Constructor for the BuildScript
    #
    #   NOTE: You can use this to add in any additional
    #   build logic you want at this stage, or add in
    #   any additional objects you want.
    #
    #   return -> BuildScript (instance of this class)
    #
    #   @param targets A dict object like the example at the top
    def __init__(self, targets):

        ### You can modify stuff in here
        self._targets = targets.copy()
        self._total_time = 0

        # Iterate through map to create a tree for the build
        self._tree = {}
        # Note the level where we are in the tree
        level = 0
        # Keep a running list of what has been 'built' already
        built = []
        while len(targets) > 0:
            # Create blank tree level
            self._tree[level] = []

            # Find build targets for current tree level
            for name in targets:
                if level == 0:
                    if targets[name]['dependencies'] == None:
                        self._tree[level].append(name)
                        built.append(name)
                else:
                    # Find whose dependencies have been previously satisfied
                    if set(targets[name]['dependencies']).issubset(built):
                        self._tree[level].append(name)
                        built.append(name)

            # Remove those targets from the possible targets
            for name in self._tree[level]:
                targets.pop(name)

            # Update the level in the tree
            level += 1
            # End while loop


    # Function to perform the simulated build
    #
    #   Your function should call the 'build' function from the
    #   'BuildTarget' helper class to specify what target is being
    #   currently built and how long it takes to build. You should
    #   also be able to supress this output.
    #
    #   return -> None
    #
    #   @param quite When set to 'True', should suppress output from 'build'
    #   @param use_parallel (Optional) When set to 'True', should build as though done in parallel
    #
    #   NOTE: You are not expected to actually implement any parallelism here, just
    #   a calculation of what the expected build time would be if using maximum/perfect
    #   parallelism when building the targets that can be built all at same time. This
    #   is a theoretical excerise here.
    def build(self, quite=False, use_parallel = False):
        ### You can modify stuff in here

        # Use the tree object to determine the actual build
        self._total_time = 0
        max_time = 0
        for level in self._tree:
            max_time = 0
            for name in self._tree[level]:
                obj = BuildTarget(name, self._targets[name]['build_time'])
                if not quite:
                    obj.build()
                if use_parallel == False:
                    self._total_time = self._total_time + self._targets[name]['build_time']
                else:
                    if self._targets[name]['build_time'] > max_time:
                        max_time = self._targets[name]['build_time']
            if use_parallel:
                self._total_time = self._total_time + max_time

    # Function to calculate sequential build time
    #
    #   return -> time (int)
    #
    #   @param quite When set to 'True', should suppress output from 'build'
    def calculate_build_time(self, quite=True):
        ### You can modify stuff in here

        # Call build function and return total time
        self.build(quite)
        return self._total_time

    # If you could build multiple targets at the same time,
    # what would the build time be?
    #
    #   return -> time (int)
    #
    #   @param quite When set to 'True', should suppress output from 'build'
    def calculate_build_time_parallel(self, quite=True):
        ### You can modify stuff in here

        # Call build function and return total time
        self.build(quite,use_parallel=True)
        return self._total_time


if __name__=="__main__":
    test = BuildScript(targets)

    ## Test 1: Make sure you have 'roughly' the correct build order
    test.build()

    # should print something to the order of (exact order can vary)
    # Built Target: 'c' Took: 5s
    # Built Target: 'e' Took: 9s
    # Built Target: 'f' Took: 8s
    # Built Target: 'g' Took: 7s
    # Built Target: 'd' Took: 2s
    # Built Target: 'b' Took: 1s
    # Built Target: 'a' Took: 10s



    ## Test 2: Calculate sequential build time
    time = test.calculate_build_time()
    if time == 42:
        print("Test 2: PASSED")
    else:
        print("Test 2: FAILED")
        print("Your Build Time: " + str(time))
        print("Expected Build Time: " + str(42))


    ## Test 3: Calculate parallel build time
    time = test.calculate_build_time_parallel()
    if time == 22:
        print("Test 3: PASSED")
    else:
        print("Test 3: FAILED")
        print("Your Build Time: " + str(time))
        print("Expected Build Time: " + str(42))
