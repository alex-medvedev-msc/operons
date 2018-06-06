from gensim.models import Doc2Vec


def main():
    doc2vec = Doc2Vec.load('models/default_tagged.doc2vec')
    similar = doc2vec.docvecs.most_similar('WP_000822720.1')
    for doc in similar:
        print(doc)
    print(doc2vec)

if __name__ == '__main__':
    main()
