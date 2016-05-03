#include "dkl_args.h"

void __args_set_field(size_t argc, char **argv, size_t n_keywords, char **keywords, bool required, void *field_ptr, char type) {
    bool found = false;
    dkl_errno = 0;
    char *err = NULL;
    char *keyword_name = keywords[n_keywords - 1];
    for(size_t i = 0; i < argc; i ++) {
        for(size_t j = 0; j < n_keywords; j ++) {
            if (streq(argv[i], keywords[j])) {
                found = true;
                /* Check boolean */
                if (type == 'b') {
                    bool *field = (bool *)field_ptr;
                    *field = true;
                } else {
                    /* Check that there are enough arguments */
                    if (argc <= i + 1) {
                        error_print("Argument value for \"%s\" is missing", keyword_name);
                        exit(DKL_ARG_PARSE_ERROR);
                    }
                    if(type == 'd') {
                        double *field = (double *)field_ptr;
                        *field = strtod(argv[i+1], &err);
                        if (argv[i+1] == err) {
                            error_print("Could not parse the value of \"%s\"", keyword_name);
                            exit(DKL_ARG_PARSE_ERROR);
                        }
                    } else if (type == 'i') {
                        int64_t *field = (int64_t *)field_ptr;
                        *field = strtol(argv[i+1], &err, 10);
                        if (argv[i+1] == err) {
                            error_print("Could not parse the value of \"%s\"", keyword_name);
                            exit(DKL_ARG_PARSE_ERROR);
                        }
                    } else if (type == 's') {
                        char **field = (char **)field_ptr;
                        *field = argv[i+1];
                    } else {
                        error_print("Unknown type in __args_set_field");
                        exit(DKL_RUNTIME_ERROR);
                    }
                }
            }
        }
    }
    if (required && !found) {
        error_print("Could not find the required option \"%s\"", keyword_name);
        exit(DKL_ARG_PARSE_ERROR);
    }
    if (!required && !found) {
        dkl_errno = DKL_OPTION_NOT_FOUND;
    }
}

bool dkl_args_parse_flag(size_t argc, char **argv, bool required, ...) {
    char **keywords = dkl_alloc(1024, char *);
    char *word;
    size_t n_keywords = 0;
    va_list words;
    va_start(words, required);
    word = va_arg(words, char *);
    while (word != NULL) {
        size_t word_size = strlen(word) + 1;
        keywords[n_keywords] = dkl_alloc(word_size, char);
        strcpy(keywords[n_keywords], word);
        n_keywords++;
        word = va_arg(words, char *);
    }
    va_end(words);
    keywords = dkl_realloc(keywords, n_keywords, char *);

    bool field = false;
    __args_set_field(argc, argv, n_keywords, keywords, required, (void *)&field, 'b');

    for (size_t i = 0; i < n_keywords; i++) {
        dkl_del(keywords[i]);
    }
    dkl_del(keywords);

    return field;
}

double dkl_args_parse_double(size_t argc, char **argv, bool required, ...) {
    char **keywords = dkl_alloc(1024, char *);
    char *word;
    size_t n_keywords = 0;
    va_list words;
    va_start(words, required);
    word = va_arg(words, char *);
    while (word != NULL) {
        size_t word_size = strlen(word) + 1;
        keywords[n_keywords] = dkl_alloc(word_size, char);
        strcpy(keywords[n_keywords], word);
        n_keywords++;
        word = va_arg(words, char *);
    }
    va_end(words);
    keywords = dkl_realloc(keywords, n_keywords, char *);

    double field = 0;
    __args_set_field(argc, argv, n_keywords, keywords, required, (void *)&field, 'd');

    for (size_t i = 0; i < n_keywords; i++) {
        dkl_del(keywords[i]);
    }
    dkl_del(keywords);

    return field;
}

int64_t dkl_args_parse_int(size_t argc, char **argv, bool required, ...) {
    char **keywords = dkl_alloc(1024, char *);
    char *word;
    size_t n_keywords = 0;
    va_list words;
    va_start(words, required);
    word = va_arg(words, char *);
    while (word != NULL) {
        size_t word_size = strlen(word) + 1;
        keywords[n_keywords] = dkl_alloc(word_size, char);
        strcpy(keywords[n_keywords], word);
        n_keywords++;
        word = va_arg(words, char *);
    }
    va_end(words);
    keywords = dkl_realloc(keywords, n_keywords, char *);

    int64_t field = 0;
    __args_set_field(argc, argv, n_keywords, keywords, required, (void *)&field, 'i');

    for (size_t i = 0; i < n_keywords; i++) {
        dkl_del(keywords[i]);
    }
    dkl_del(keywords);

    return field;
}

char *dkl_args_parse_string(size_t argc, char **argv, bool required, ...) {
    char **keywords = dkl_alloc(1024, char *);
    char *word;
    size_t n_keywords = 0;
    va_list words;
    va_start(words, required);
    word = va_arg(words, char *);
    while (word != NULL) {
        size_t word_size = strlen(word) + 1;
        keywords[n_keywords] = dkl_alloc(word_size, char);
        strcpy(keywords[n_keywords], word);
        n_keywords++;
        word = va_arg(words, char *);
    }
    va_end(words);
    keywords = dkl_realloc(keywords, n_keywords, char *);

    char *field = dkl_new(char *);
    field = NULL;
    __args_set_field(argc, argv, n_keywords, keywords, required, (void *)&field, 's');

    for (size_t i = 0; i < n_keywords; i++) {
        dkl_del(keywords[i]);
    }
    dkl_del(keywords);

    return field;
}
