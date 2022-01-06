from dataforest.hooks import hook


@hook
def hook_example_no_args(dp):
    print("hook_example_no_args")


@hook("requires")
def hook_example_args(dp):
    print("hook_example_args")
