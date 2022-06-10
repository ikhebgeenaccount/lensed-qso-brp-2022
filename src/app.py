import configparser


class App:

    __conf = None

    @staticmethod
    def config():
        if App.__conf is None:  # Read only once, lazy.
            App.__conf = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation(),
                                                   converters={'list': lambda s: [st.strip() for st in s.split(',')]})
            App.__conf.optionxform = str  # To make sure options stay capitalized as they are, default makes them all lower case
            App.__conf.read('settings.ini')
        return App.__conf
