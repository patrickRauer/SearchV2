import os
from datetime import datetime

TIME_FORMAT = '%y-%m-%d %H:%M:%S'


class LogItem:
    time = None
    action = None
    comment = None

    def __init__(self, action, comment, time=None):
        """
        Log-item for the Log-class. It contains the needed information
        # TODO: small update of the description

        :param action: The performed action
        :type action: str
        :param comment: An optional comment for the action
        :type comment: str
        :param time:
            The time of the action, default is None which means
            that the current time will be used.
        :type time: datetime.datetime, None
        """
        if time is None:
            time = datetime.now()
        self.time = time
        self.action = action
        self.comment = comment

    def __str__(self):
        """
        Converts the LogItem to a string in a standardized format

        :return: The data of the LogItem as a string
        :type: str
        """
        return '{}\t{}\t{}\n'.format(self.time.strftime(TIME_FORMAT),
                                     self.action,
                                     self.comment)

    @staticmethod
    def from_str(log_str):
        """
        Converts a log-item string to a LogItem object

        :param log_str: The string with log item data
        :type log_str: str
        :return: A LogItem with the data of the string
        :rtype: Search_V2.core.log.LogItem
        """

        # prepare the input string
        log_str = log_str.strip('\n')
        log_str = log_str.split('\t')

        # set the data from the string to the LogItem variables
        log_item = LogItem(
            log_str[1],
            log_str[-1],
            datetime.strptime(log_str[0],
                              TIME_FORMAT)
        )
        return log_item


class Log:
    items = []

    def __init__(self):
        """
        This class provides a log system to track the process of the search.
        """
        pass

    def add_log(self, action, comment=''):
        """
        Adds a new log item to the log.

        :param action: The performed action
        :type action: str
        :param comment: An optional comment of the action
        :type comment: str
        :return:
        """
        self.items.append(LogItem(action, comment))

    def write(self, path, overwrite=False):
        """
        Writes the log data to a file.

        :param path: Path to the file
        :type path: str
        :param overwrite:
            True if an existing file should be overwritten, else False. Default is false and if a file
            at the place exists, an error will raise.
        :return:
        """
        if os.path.exists(path) and not overwrite:
            raise FileExistsError('File already exists!')

        # check if the directory doesn't exists and created if needed
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        # write all the log-items to the file
        with open(path, 'w') as f:
            for i in self.items:
                f.write(str(i))

    def __str__(self):
        out = ''

        for i in self.items:
            out += str(i)
        return out

    @staticmethod
    def load(path):
        """
        Load an existing log-file

        :param path: Path to the file
        :type path: str
        :return: The log-file as a Log-object
        :rtype: Search_V2.core.log.Log
        """
        log = Log()
        with open(path) as f:
            for line in f:
                log.items.append(LogItem.from_str(line))
        return log
