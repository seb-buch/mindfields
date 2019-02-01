from collections import OrderedDict
import jsonlines
import re
import sys
import os
from typing import Union
import argparse


class CorpusEntry(object):
    def __init__(self,
                 text: str,
                 meta: OrderedDict = None):
        self.text = text
        if meta is None:
            meta = OrderedDict()
        self.meta = OrderedDict(meta)

    def add_meta(self, key: str, value: str):
        self.meta[key] = value

    def as_dict(self):

        return {
            "text": self.text,
            "meta": self.meta
        }


def is_useful_for_corpus(text: Union[str, bytes]) -> bool:
    required_words = ["anti", "inhibit"]

    if type(text) == bytes:
        text = text.decode("utf-8")

    for word in required_words:
        if word in text:
            return True

    return False


def already_in_corpus(text: str, corpus: [CorpusEntry]) -> bool:
    for entry in corpus:
        if entry.text == text:
            return True

    return False


def get_uniprot_entry_length(xml_entry: bytes) -> int:
    start_pos = xml_entry.find(b"<sequence")

    if start_pos < 0:
        return -1

    start_pos = xml_entry.find(b"length=", start_pos)

    if start_pos < 0:
        return -1

    attribute = xml_entry[start_pos:].split()[0]
    name, value = attribute.split(b"=")

    return int(value.strip(b'"'))


def populate_uniprot_corpus(corpus: list = None,
                            sprot_file: str = "data/uniprot_sprot.xml",
                            max_added: int = 1e6,
                            max_length: int = 100) -> [CorpusEntry]:

    if corpus is None:
        corpus = list()

    corpus_size = len(corpus)
    size = os.path.getsize(sprot_file)
    buffer_size = int(10e6)

    start_tag = b"<entry "
    end_tag = b"</entry>"

    id_pattern = re.compile(b"<accession[^>]*>(?P<uniprot_id>[^<]+)</accession>")

    comment_text_pattern = re.compile(b"<text[^>]*>(?P<text>[^<]+)</text>")

    article_title_pattern = re.compile(b"<title[^>]*>(?P<text>[^<]+)</title>")

    nused = 0
    ntotal = 0
    nuniprot_entries = 0

    try:
        with open(sprot_file, "rb") as fp:
            current_position = fp.tell()

            while current_position < size:
                print("Reading '{}' (buffer size: {:.1f} MB): {:5.1f}% Done ({} uniprot entries read ->"
                      "{}/{} corpus entries used)\r".format(
                    sprot_file,
                    buffer_size * 1e-6,
                    current_position / size * 100,
                    nuniprot_entries,
                    nused,
                    ntotal
                ), end="")
                if nused >= max_added:
                    print("\nLimit reached for the number of entries added to corpus", end="")
                    break
                sys.stdout.flush()
                current_position = fp.tell()

                buffer = fp.read(buffer_size)

                begin_pos = buffer.find(start_tag)
                counter = 0
                while begin_pos > -1:
                    end_pos = buffer.find(end_tag, begin_pos)
                    if end_pos > 0:  # A full entry was found
                        nuniprot_entries += 1
                        counter += 1

                        xml_entry = buffer[begin_pos:end_pos+len(end_tag)]

                        if get_uniprot_entry_length(xml_entry) < max_length:
                            # Only consider entries that correspond to peptides or small proteins

                            meta = OrderedDict()
                            meta["source"] = "Uniprot"

                            # Add meta data for uniprotid
                            for match in id_pattern.finditer(xml_entry):
                                meta["ID"] = match.group("uniprot_id").strip().decode("utf-8")

                            for match in comment_text_pattern.finditer(xml_entry):

                                sentences = match.group("text").decode("utf-8").split(". ")

                                for text in sentences:
                                    ntotal += 1

                                    if is_useful_for_corpus(text) and not already_in_corpus(text, corpus):
                                        nused += 1
                                        if nused > max_added:
                                            break

                                        entry = CorpusEntry(text, meta=meta)
                                        entry.add_meta("type", "comment")

                                        corpus.append(entry)

                            for match in article_title_pattern.finditer(xml_entry):
                                sentences = match.group("text").decode("utf-8").split(". ")

                                for text in sentences:
                                    ntotal += 1

                                    if is_useful_for_corpus(text) and not already_in_corpus(text, corpus):
                                        nused += 1
                                        if nused > max_added:
                                            break

                                        entry = CorpusEntry(text, meta=meta)
                                        entry.add_meta("type", "article")

                                        corpus.append(entry)

                        # We set the file pointer to the end of xml entry
                        fp.seek(current_position + end_pos + len(end_tag))

                    elif counter == 0:
                        raise ValueError("Buffer size is too small to get a full entry!")
                    else:
                        fp.seek(current_position + begin_pos)
                        break

                    if nused >= max_added:
                        break
                    begin_pos = buffer.find(start_tag, begin_pos+1)

    except KeyboardInterrupt:
        print("\nInterruption requested!")
    else:
        print("")

    print("INFO: {} entries added to corpus".format(len(corpus) - corpus_size))
    return corpus


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--max-size", help="Corpus maximum size", type=int, default=1000)

    parser.add_argument("--max-length", help="Maximum sequence length", type=int, default=100)

    parser.add_argument("--output", "-o", help="Name for the .jsonl file used to store the corpus",
                        default="../raw_corpus.jsonl")

    args = parser.parse_args()

    corpus = populate_uniprot_corpus(max_added=args.max_size, max_length=args.max_length)

    with jsonlines.open(args.output, mode='w', compact=True) as writer:
        for entry in corpus:
            writer.write(entry.as_dict())

    print("INFO: Corpus saved to '{}'".format(args.output))
