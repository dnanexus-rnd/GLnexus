// This is an example adapted from the tutorial
// (https://capnproto.org/cxx.html).

// C files handling file-descriptors
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "tut_proto.capnp.h"
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <iostream>

using namespace GLnexus::serialization;

void writeAddressBook(int fd) {
    ::capnp::MallocMessageBuilder message;

    AddressBook::Builder addressBook = message.initRoot<AddressBook>();
    ::capnp::List<Person>::Builder people = addressBook.initPeople(2);

    Person::Builder alice = people[0];
    alice.setId(123);
    alice.setName("Alice");
    alice.setEmail("alice@example.com");
    // Type shown for explanation purposes; normally you'd use auto.
    ::capnp::List<Person::PhoneNumber>::Builder alicePhones =
          alice.initPhones(1);
    alicePhones[0].setNumber("555-1212");
    alicePhones[0].setType(Person::PhoneNumber::Type::MOBILE);
    alice.getEmployment().setSchool("MIT");

    Person::Builder bob = people[1];
    bob.setId(456);
    bob.setName("Bob");
    bob.setEmail("bob@example.com");
    auto bobPhones = bob.initPhones(2);
    bobPhones[0].setNumber("555-4567");
    bobPhones[0].setType(Person::PhoneNumber::Type::HOME);
    bobPhones[1].setNumber("555-7654");
    bobPhones[1].setType(Person::PhoneNumber::Type::WORK);
    bob.getEmployment().setUnemployed();

    writePackedMessageToFd(fd, message);
}

void printAddressBook(int fd) {
    ::capnp::PackedFdMessageReader message(fd);

    AddressBook::Reader addressBook = message.getRoot<AddressBook>();

    for (Person::Reader person : addressBook.getPeople()) {
        std::cout << person.getName().cStr() << ": "
                  << person.getEmail().cStr() << std::endl;
        for (Person::PhoneNumber::Reader phone: person.getPhones()) {
            const char* typeName = "UNKNOWN";
            switch (phone.getType()) {
            case Person::PhoneNumber::Type::MOBILE: typeName = "mobile"; break;
            case Person::PhoneNumber::Type::HOME: typeName = "home"; break;
            case Person::PhoneNumber::Type::WORK: typeName = "work"; break;
            }
            std::cout << "  " << typeName << " phone: "
                      << phone.getNumber().cStr() << std::endl;
        }
        Person::Employment::Reader employment = person.getEmployment();
        switch (employment.which()) {
        case Person::Employment::UNEMPLOYED:
            std::cout << "  unemployed" << std::endl;
            break;
        case Person::Employment::EMPLOYER:
            std::cout << "  employer: "
                      << employment.getEmployer().cStr() << std::endl;
            break;
        case Person::Employment::SCHOOL:
            std::cout << "  student at: "
                      << employment.getSchool().cStr() << std::endl;
            break;
        case Person::Employment::SELF_EMPLOYED:
            std::cout << "  self-employed" << std::endl;
            break;
        }
    }
}

int main(int argc, char* argv[]) {
    char *filename = (char*) "xyz";

    {
        int fd = open(filename, O_WRONLY | O_TRUNC | O_CREAT, 0644);
        if (fd < 0) {
            fprintf(stderr, "Could not open file %s for writing\n", filename);
            return 1;
        }
        writeAddressBook(fd);
        close(fd);
    }

    {
        int fd = open(filename, O_RDONLY);
        if (fd < 0) {
            fprintf(stderr, "Could not open file %s for reading\n", filename);
            return 1;
        }
        printAddressBook(fd);
        close(fd);
    }

    return 0;
}
