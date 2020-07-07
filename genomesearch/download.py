from genomesearch import *

def _download():
    click.echo("#### INPUT PARAMETERS ####")
    try:
        num_markers = int(input("How many markers do you want to use? (This can be any number between 1 and 400, default is 150)") or "150")
        if num_markers > 400 or num_markers < 1:
            raise Exception('wrong_markers')
    except:
        print("ERROR!")
        print("Please input a number between 1 and 400, the default is 150.")
        num_markers = int(input("How many markers do you want to use?") or "150")

    print(num_markers)
    click.echo("####################")