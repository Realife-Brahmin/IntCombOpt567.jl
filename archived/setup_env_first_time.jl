using PkgTemplates

tpl = Template(
    user="Realife-Brahmin",              # or any user/organization name you prefer
    dir = pwd(), # will create the folder IntCombOpt567 in the current directory which basically is the home of the package
    authors=["Aryan Ritwajeet Jha <aryan.r.jha@gmail.com>"], # optional
    license="MIT",                          # or another license name
    julia_compat="1.6",                     # or whichever versions you want to support
    plugins=[
        # You can add or remove plugins as you see fit
    ]
)

generate("IntCombOpt567", tpl)